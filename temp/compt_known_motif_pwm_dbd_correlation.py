#!/usr/bin/python
import os
import sys
import numpy
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

# species = "yeast"; species_full = "Saccharomyces_cerevisiae"
species = "fly"; species_full = "Drosophila_melanogaster"
min_score_dbd = 0
itv_score_dbd = 10
# score_type = "median"
score_type = "average"
parse_method = "e_value"
# parse_method = "-log10_fc"
fit_function = "sigmoid"
# fit_function = "exponential"

def main():
	# dir_proj = "/home/mblab/ykang/proj_db_infer_pipe/"
	dir_proj = "/Users/KANG/cgscluster/proj_db_infer_pipe/"
	fn_motif_list = dir_proj + "resources/cisbp_1.01/cisbp_motifs_" + species_full + ".txt"
	dir_dbd = dir_proj + "output/" + species + "_cisbp_known_motif_dbd_pid/"
	dir_pwm = dir_proj + "output/" + species + "_cisbp_known_motif_pwm_similarity/"

	motifs = numpy.loadtxt(fn_motif_list, dtype=str)

	# parse data
	scores_dbd = []
	scores_pwm = []
	for motif in motifs:
		fn_dbd = dir_dbd + motif
		fn_pwm = dir_pwm + motif + ".tomtom"
		if (os.path.isfile(fn_dbd) and os.path.getsize(fn_dbd) > 0) and (os.path.isfile(fn_pwm) and os.path.getsize(fn_pwm) > 0):
			dict_dbd = parse_dict(fn_dbd, True, min_score_dbd, None)
			dict_pwm = parse_dict(fn_pwm, False, None, parse_method)
			if bool(dict_dbd) and bool(dict_pwm):
				targets = numpy.intersect1d(dict_dbd.keys(), dict_pwm.keys())
				targets_used = []
				for target in targets:
					if not target in targets_used:
						scores_dbd.append(dict_dbd[target])
						scores_pwm.append(dict_pwm[target])
						targets_used.append(target)

	# bin scores
	[bins_scores_dbd, med_scores_pwm, avg_scores_pwm] = compt_bins(scores_dbd, scores_pwm, min_score_dbd, itv_score_dbd)	
	if score_type == 'median':
		y = med_scores_pwm
		label_scatter = "median scores"
	elif score_type == 'average':
		y = avg_scores_pwm
		label_scatter = "average scores"
	else:
		sys.exit("Median or average score not specified")

	# logistic fit
	p_guess=(numpy.median(bins_scores_dbd), numpy.median(y), 1, .25)
	if fit_function == 'sigmoid':
		p, cov, infodict, mesg, ier = leastsq(sigmoid_residuals, p_guess, args=(bins_scores_dbd, y), full_output=1)
	elif fit_function == 'exponential':
		p, cov, infodict, mesg, ier = leastsq(exponential_residuals, p_guess, args=(bins_scores_dbd, y), full_output=1)
	else:
		sys.exit("Fit function not specified")

	# plot results
	plt.figure()
	plt.scatter(scores_dbd, scores_pwm, alpha=.5)
	plt.scatter(bins_scores_dbd, y, alpha=.5, s=100, c='g', label=label_scatter)
	xp = numpy.linspace(min_score_dbd, 100, 100)
	if fit_function == 'sigmoid':
		pxp = sigmoid(p, xp)
		label_fit = "sigmoid fit, c=" + str(p[2]) + ", k=" + str(p[3])
	elif fit_function == 'exponential':
		pxp = exponential(p, xp)
		label_fit = "exponential fit, c=" + str(p[2]) + ", k=" + str(p[3])
	else:
		sys.exit("Fit function not specified")
	plt.plot(xp, pxp, 'r-', label=label_fit)
	plt.legend(loc="lower right")
	plt.xlim([0, 100])
	plt.ylim([0, 1])
	plt.title("DBD vs PWM Similarity, Species: " + species)
	plt.xlabel("DBD Percent Identity")
	plt.ylabel("PWM Similarity (1 - E_value)")
	plt.show()

def parse_dict(fn, use_cutoff, min_score_dbd, parse_method):
	lines = open(fn, 'r').readlines()
	dict = {}
	if len(lines) > 1:
		if use_cutoff:
			for line in lines:
				target = line.split()[0].split(":")[0]
				score = float(line.split()[1])
				if score < min_score_dbd:
					break
				dict[target] = score
		else:
			if parse_method == 'e_value':
				for i in range(1, len(lines)):
					target = lines[i].split()[0].split(":")[0]
					score = float(lines[i].split()[1])
					dict[target] = 1 - score
			elif parse_method == '-log10_fc':
				baseline = -numpy.log10(float(lines[0].split()[1]))
				for i in range(1, len(lines)):
					target = lines[i].split()[0].split(":")[0]
					score = -numpy.log10(float(lines[i].split()[1]))
					dict[target] = score/baseline
			else:
				sys.exit("PWM parse method not specified")
	return dict

def compt_bins(scores_dbd, scores_pwm, min_score_dbd, itv_score_dbd):
	scores_dbd = numpy.array(scores_dbd)
	scores_pwm = numpy.array(scores_pwm)
	num_bins = int((100 - min_score_dbd)/itv_score_dbd)
	bins = numpy.arange(min_score_dbd + itv_score_dbd, 100 + itv_score_dbd, itv_score_dbd)
	med_scores_pwm = numpy.zeros(num_bins)
	avg_scores_pwm = numpy.zeros(num_bins)
	for i in range(num_bins):
		temp_indices = numpy.where((scores_dbd > bins[i]-itv_score_dbd) & (scores_dbd <= bins[i]))[0]
		temp_scores = scores_pwm[temp_indices]
		med_scores_pwm[i] = numpy.median(temp_scores)
		avg_scores_pwm[i] = numpy.average(temp_scores)
	return numpy.array([bins, med_scores_pwm, avg_scores_pwm])

def sigmoid(p,x):
    x0,y0,c,k = p
    y = c / (1 + numpy.exp(-k * (x-x0))) + y0
    return y

def sigmoid_residuals(p,x,y):
	return y - sigmoid(p, x)

def exponential(p,x):
	x0,y0,c,k = p
	y = c * numpy.exp(k * (x-x0)) + y0
	return y

def exponential_residuals(p,x,y):
	return y - exponential(p, x)

if __name__ == "__main__":
    main()