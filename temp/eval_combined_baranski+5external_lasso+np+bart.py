#!/usr/bin/python

"""
Combine mutliple tf combination of the edge scores at difference dbd cutoff.

eval_method: cumulative, binned
"""

import sys
import argparse
import glob
import os.path
import numpy
import matplotlib.pyplot as plt

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine chip and pwm evaluations")
    parser.add_argument('-e', '-eval_method', dest='eval_method', type=str, default='binned')
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    # original np and bartNp evaluations
    # dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_combined/'
    # dir_sub_gmean = 'fly_network_baranski_flyat_grav_postfast_sexan_toxic_geometric_mean/analysis_flynet/'
    # dir_sub_quantile = 'fly_network_baranski_flyat_grav_postfast_sexan_toxic_quantile_combine/analysis_flynet/'
    # fns = []
    # fns.append(dir_sub_gmean + 'analysis_chip_support.lasso_only.txt')
    # fns.append(dir_sub_gmean + 'analysis_pwm_support.lasso_only.txt')
    # fns.append(dir_sub_quantile + 'analysis_chip_support.lasso_only.txt')
    # fns.append(dir_sub_quantile + 'analysis_pwm_support.lasso_only.txt')
    # fns.append(dir_sub_gmean + 'analysis_chip_support.lasso_np.txt')
    # fns.append(dir_sub_gmean + 'analysis_pwm_support.lasso_np.txt')
    # fns.append(dir_sub_quantile + 'analysis_chip_support.lasso_np.txt')
    # fns.append(dir_sub_quantile + 'analysis_pwm_support.lasso_np.txt')
    # fns.append(dir_sub_gmean + 'analysis_chip_support.lasso_np_bart.txt')
    # fns.append(dir_sub_gmean + 'analysis_pwm_support.lasso_np_bart.txt')
    # fns.append(dir_sub_quantile + 'analysis_chip_support.lasso_np_bart.txt')
    # fns.append(dir_sub_quantile + 'analysis_pwm_support.lasso_np_bart.txt')
    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_combined/'
    dir_sub_gmean = 'fly_network_baranski_flyat_grav_postfast_sexan_toxic_geometric_mean/analysis_new/'
    dir_sub_quantile = 'fly_network_baranski_flyat_grav_postfast_sexan_toxic_quantile_combine/analysis_new/'
    fns = []
    # eval_range = ['full', '40k']
    eval_range = ['genome_wide', '100k']
    fns.append(dir_sub_gmean + 'analysis_chip_support.' + eval_range[0] + '.lasso_np.txt')
    fns.append(dir_sub_gmean + 'analysis_pwm_support.' + eval_range[1] + '.lasso_np.txt')
    fns.append(dir_sub_quantile + 'analysis_chip_support.' + eval_range[0] + '.lasso_np.txt')
    fns.append(dir_sub_quantile + 'analysis_pwm_support.' + eval_range[1] + '.lasso_np.txt')
    fns.append(dir_sub_gmean + 'analysis_chip_support.' + eval_range[0] + '.lasso_np_bart.txt')
    fns.append(dir_sub_gmean + 'analysis_pwm_support.' + eval_range[1] + '.lasso_np_bart.txt')
    fns.append(dir_sub_quantile + 'analysis_chip_support.' + eval_range[0] + '.lasso_np_bart.txt')
    fns.append(dir_sub_quantile + 'analysis_pwm_support.' + eval_range[1] + '.lasso_np_bart.txt')

    # figure setup
    # colors = ['k:', 'r', 'r--', 'g', 'g--', 'b', 'b--']
    # labels = []
    # labels.append('chance')
    # labels.append('lasso: geo mean')
    # labels.append('lasso: quant comb')
    # labels.append('lasso+np: geo mean')
    # labels.append('lasso+np: quant comb')
    # labels.append('lasso+np+bart: geo mean')
    # labels.append('lasso+np+bart: quant comb')
    # x_ticks = ['2k', '4k', '6k', '8k', '10k', '12k', '14k', '16k', '18k', '20k']
    colors = ['k:', 'r', 'r--', 'b', 'b--']
    labels = []
    labels.append('chance')
    labels.append('lasso+np: geo mean')
    labels.append('lasso+np: quant comb')
    labels.append('lasso+np+bart: geo mean')
    labels.append('lasso+np+bart: quant comb')
    x_ticks = ['2k', '4k', '6k', '8k', '10k', '12k', '14k', '16k', '18k', '20k']

    # compute chip and pwm support ratios of each method
    eval_chip = numpy.zeros([len(fns)/2+1, 10])
    eval_pwm = numpy.zeros([len(fns)/2+1, 10])

    for i in range(len(fns)/2):
        chip_support = numpy.loadtxt(dir_network + fns[i*2])
        pwm_support = numpy.loadtxt(dir_network + fns[i*2+1]) 
        if i == 0:
            eval_chip[0,:] = chip_support[0,:]
            eval_pwm[0,:] = pwm_support[0,:]
        eval_chip[i+1,:] = chip_support[1,:]
        eval_pwm[i+1,:] = pwm_support[1,:]

    if parsed.eval_method == 'cumulative':
        temp_eval_chip = eval_chip
        temp_eval_pwm = eval_pwm
        for i in range(10):
            eval_chip[:,i] = numpy.sum(temp_eval_chip[:,0:(i+1)], axis=1)/(i+1)
            eval_pwm[:,i] = numpy.sum(temp_eval_pwm[:,0:(i+1)], axis=1)/(i+1)

    # plot figures
    plt.figure(num=None, figsize=(15,8), dpi=80)
    plt.subplot(1,2,1)
    for i in range(len(eval_chip)):
        plt.plot(eval_chip[i], colors[i], label=labels[i])
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by ChIP')
    plt.xlim(-1, len(eval_chip[0]))
    plt.ylim(0, 0.45)
    plt.legend(loc="upper right")

    plt.subplot(1,2,2)
    for i in range(len(eval_pwm)):
        plt.plot(eval_pwm[i], colors[i], label=labels[i])
    plt.xticks(range(len(eval_pwm[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by PWM')
    plt.xlim(-1, len(eval_pwm[0]))
    plt.ylim(0, 0.15)
    plt.legend(loc="upper right")

    plt.show()

def parse_ratio(fn, method):
    lines = open(fn, "r").readlines()
    chip = [0] * (len(lines)/3)
    pwm = [0] * (len(lines)/3)
    if method == "cumulative":
        for i in range(len(lines)/3):
            line = lines[i*3].split()
            chip[i] = float(line[5])/float(line[2])
            pwm[i] = float(line[4])/float(line[2])
    elif method == "binned":
        for i in range(len(lines)/3):
            line = lines[i*3].split()
            if i == 0:
                chip[i] = float(line[5])/float(line[2])
                pwm[i] = float(line[4])/float(line[2])
            else:
                chip[i] = (float(line[5]) - float(prevline[5]))/(float(line[2]) - float(prevline[2]))
                pwm[i] = (float(line[4]) - float(prevline[4]))/(float(line[2]) - float(prevline[2]))
            prevline = line    
    return  [chip, pwm]

if __name__ == "__main__":
    main(sys.argv)
