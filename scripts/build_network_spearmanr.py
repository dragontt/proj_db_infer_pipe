#!/usr/bin/python
import sys
import argparse
import numpy
from scipy.stats import spearmanr

def parse_args(argv):
	parser = argparse.ArgumentParser(description="Build a network using Spearman correlation of expression data")
	parser.add_argument('-r', '--expr_regulator', dest='expr_regulator', type=str)
	parser.add_argument('-g', '--expr_gene', dest='expr_gene', type=str)
	parser.add_argument('-a', '--allowed', dest='allowed', type=str)
	parser.add_argument('-o', '--output_adjmtr', dest='output_adjmtr', type=str)
	parsed = parser.parse_args(argv[1:])
	return parsed

def main(argv):
	parsed = parse_args(argv)

	# parse expression data
	expr_r = numpy.loadtxt(parsed.expr_regulator)
	expr_g = numpy.loadtxt(parsed.expr_gene)
	allowed = numpy.loadtxt(parsed.allowed)
	(n,m) = allowed.shape

	# disallow the genes with all zero expressions
	for i in range(n):
		if sum(expr_r[i]) == 0:
			allowed[i,:] = 0
	for i in range(m):
		if sum(expr_g[i]) == 0:
			allowed[:,i] = 0

	# build the network using spearman correlation
	adjmtr = numpy.zeros((n,m))
	for i in range(n):
		for j in range(m):
			if allowed[i,j] == 1:
				corr = spearmanr(expr_r[i], expr_g[j])
				adjmtr[i,j] = corr[0]

	# write data
	numpy.savetxt(parsed.output_adjmtr, adjmtr, fmt='%.5f', delimiter='\t', newline='\n')

if __name__ == "__main__":
	main(sys.argv)
