#!/usr/bin/python
""" 
Compute the fold change of expression data. 
Methods:
- wt
Use this method if the expression data contain gene perturbation (overexpression or deletion). Duplicated WT samples are averaged. Pseudocount is set at 0.5 percentile of all nonzero expression data from WT and non-WT.
- mean
Use this method if the expression data do NOT contain gene perturbation. The baseline for each gene is the geometric mean of the expression profiles in all conditions. 
"""

import sys
import argparse
import numpy
from scipy.stats import gmean

def parse_args(argv):
	parser = argparse.ArgumentParser(description="Compute fold change of expression data. Use method 'wt' (defualt) if the expression data contain gene perturbation; use method 'mean' if the expression data do NOT contain gene perturbation.")
	parser.add_argument('-m', '--method', dest='method', default='wt')
	parser.add_argument('-e', '--fn_expr', dest='fn_expr')
	parser.add_argument('-c', '--fn_cond', dest='fn_cond')
	parser.add_argument('-g', '--fn_gids', dest='fn_gids')
	parser.add_argument('-o', '--fn_output', dest='fn_output')
	parsed = parser.parse_args(argv[1:])
	return parsed

def main(argv):
	parsed = parse_args(argv)

	# load data
	expr = numpy.transpose(numpy.loadtxt(parsed.fn_expr))
	cond = numpy.loadtxt(parsed.fn_cond, dtype=str, delimiter='\t')
	gids = numpy.loadtxt(parsed.fn_gids, dtype=str)

	# Method - wt
	if parsed.method == 'wt':
		# compute pseudocount
		expr_indices = numpy.transpose(numpy.nonzero(expr))
		expr_nonzero = numpy.zeros(len(expr_indices))
		for i in range(expr_indices.shape[0]):
			expr_nonzero[i] = expr[expr_indices[i,0], expr_indices[i,1]]
		pseudocount = numpy.percentile(expr_nonzero, 0.5)
		print pseudocount
		expr += pseudocount

		# expression with specified sample environmental condition
		if len(cond[0].split('-')) > 1:
			# parse indices
			wt_index_groups = []
			nonwt_indices = []
			prev_cond = ''
			wt_conds = []

			for i in range(len(cond)):
				line = cond[i].split('-')
				curr_id = line[0]
				curr_cond = '-'.join(line[1:len(line)-3]) + '-' + line[len(line)-2]
				
				if curr_id == '00000':
					if curr_cond == prev_cond:
						wt_index_groups[len(wt_index_groups)-1].append(i)
					else:
						wt_index_groups.append([i])
						wt_conds.append(curr_cond)
					prev_cond = curr_cond
				else:
					nonwt_indices.append(i)

			# compute geometric mean of wt duplicates 
			expr_fc = numpy.ones(expr.shape)
			wt_expr = {}

			for i in range(len(wt_index_groups)):
				if len(wt_index_groups[i]) == 1:
					wt_expr[wt_conds[i]] = expr[wt_index_groups[i][0]]
				else:
					temp_combined_expr = [expr[wt_index_groups[i][0]]]
					for j in range(1,len(wt_index_groups[i])):
						temp_combined_expr.append(expr[wt_index_groups[i][j]])
					temp_gmean_expr = gmean(temp_combined_expr)
					wt_expr[wt_conds[i]] = temp_gmean_expr

			# compute fold change for non-wt
			for i in nonwt_indices:
				line = cond[i].split('-')
				curr_cond = '-'.join(line[1:len(line)-3]) + '-' + line[len(line)-2]
				expr_fc[i] = expr[i] / wt_expr[curr_cond]

		# expression with NO specified sample environmental condition
		elif len(cond[0].split('-')) == 1:
			# parse indices
			wt_index = []
			cond_reps = {}
			cond_out = [None] * len(cond)
			for i in range(len(cond)):
				if cond[i].strip() == 'WT':
					wt_index.append(i)
				if not cond[i] in cond_reps.keys():
					cond_reps[cond[i]] = [0]
				else:
					cond_reps[cond[i]].append(len(cond_reps[cond[i]])+1)
				cond_out[i] = cond[i] + '_rep' + str(cond_reps[cond[i]][len(cond_reps[cond[i]])-1])

			# compute geometric mean of wt duplicates
			wt_expr = gmean(expr[wt_index,:])

			# compute fold change
			expr_fc = numpy.zeros(expr.shape)
			for i in range(len(cond)):
				expr_fc[i] = expr[i] / wt_expr

	# Method - mean
	elif parsed.method == 'mean':
		cond_out = cond
		baseline_expr = gmean(expr)
		expr_fc = numpy.zeros(expr.shape)
		for i in range(len(cond)):
			expr_fc[i] = expr[i] / baseline_expr

	# write fold change data in tsv format
	writer = open(parsed.fn_output, 'w')

	for i in range(len(gids)):
		writer.write('\"%s\"\t' % gids[i])
	writer.write('\n')

	for i in range(len(expr_fc)):
		writer.write('\"%s\"\t' % (cond_out[i]))
		for j in range(len(expr_fc[i])):
			if expr_fc[i,j] == 1:
				writer.write('1\t')
			else:
				writer.write('%.10f\t' % expr_fc[i,j])
		writer.write('\n')
	writer.close()

if __name__ == "__main__":
    main(sys.argv)
