#!/usr/bin/python

import sys
import numpy
import matplotlib.pyplot as plt

def main():
	""" yeast motifs """
	# dir_figures = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/yeast_analysis_results/'
	# figure_title = 'Motif Alignment Evaluation: E-value Cutoffs'

	# filename = 'tomtom_E_value_30.txt'
	# fns = []

	# # fire inference
	# dir_output = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/yeast_motif_inference/'
	# fns.append(dir_output + 'fire_7mers_holstege_np_bart_tf_merged_dbd50/' + filename)
	# # fire ortho inference
	# dir_output = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/yeast_network_holstege_orthologs/'
	# fns.append(dir_output + 'ortho_scer+smik+skud+sbay/' + filename)
	# fns.append(dir_output + 'ortho_scer+smik+skud+sbay+scas+sklu/' + filename)
	# # cisbp inference
	# dir_output = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/yeast_motif_inference/'
	# fns.append(dir_output + 'cisbp_dbd_cutoff_holstege_np_bart_tf_merged_dbd50/' + filename)

	""" yeast motifs """
	dir_figures = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_analysis_results/'
	figure_title = 'Motif Alignment Evaluation: E-value Cutoffs'

	filename = 'tomtom_E_value_30.txt'
	fns = []

	# fire inference
	dir_output = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_motif_inference/'
	fns.append(dir_output + 'fire_7mers_cellCycle_np_tf_merged_dbd50/' + filename)
	# fire ortho inference
	dir_output = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_orthologs/'
	fns.append(dir_output + 'motif_inference_dmel+Dsim+Dsec/' + filename)
	fns.append(dir_output + 'motif_inference_dmel+Dsim+Dsec+Dyak+Dere+Dana/' + filename)
	# cisbp inference
	dir_output = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_motif_inference/'
	fns.append(dir_output + 'cisbp_-2000_+200_fimo_dbd_cutoff_cellCycle_np_tf_merged_dbd50/' + filename)

	# compute histogram
	cutoffs = [.25,.5,.75,1,10,20,30]
	counts = []
	for i in range(len(fns)):
		counts.append(compute_histogram(fns[i], cutoffs))

	# plot bar plot
	indices = range(len(cutoffs))
	width = .2
	fig, ax = plt.subplots(figsize=(8,8), dpi=80)
	bars = []
	colors = ['r', 'g', 'c', 'b']
	legend = ('fire_motif', 'fire_motif_3_ortho_species', 'fire_motif_5_ortho_species', 'cisbp_motif')
	for i in range(len(counts)):
		bars.append(ax.bar([x + width*i for x in indices], counts[i], width, color=colors[i], alpha=.5))
	ax.set_xticks([x + width*(len(cutoffs)/2) for x in indices])
	ax.set_xticklabels(tuple([str(x) for x in cutoffs]))
	ax.set_xlabel('E-value Cutoffs')
	ax.set_ylabel('Counts of Motif Alignment E-value < Cutoff')
	ax.set_title(figure_title)
	ax.legend(legend)

	plt.savefig(dir_figures + 'motif_eval.pdf', fmt='pdf')
	plt.show()

def compute_histogram(fn, cutoffs):
	cutoffs = [0] + cutoffs
	data = numpy.loadtxt(fn, skiprows=1, usecols=[7])
	counts = []
	for i in range(1,len(cutoffs)):
		counts.append(len(numpy.where((data < cutoffs[i]) & (data >= cutoffs[i-1]))[0]))
	return counts

if __name__ == "__main__":
    main()