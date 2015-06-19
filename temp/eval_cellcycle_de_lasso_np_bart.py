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

    """ evaluate chip and pwm supports on binary gold standard """
    # # file initialization
    # dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_combined/'
    # dir_sub_quantile = 'quantile_combine_cellCycle/analysis_flynet/'
    # fns = []
    # fns.append(dir_network + dir_sub_quantile + 'analysis_chip_support.de.txt')
    # fns.append(dir_network + dir_sub_quantile + 'analysis_pwm_support.de.txt')
    # fns.append(dir_network + dir_sub_quantile + 'analysis_chip_support.lasso.txt')
    # fns.append(dir_network + dir_sub_quantile + 'analysis_pwm_support.lasso.txt')
    # fns.append(dir_network + dir_sub_quantile + 'analysis_chip_support.np.txt')
    # fns.append(dir_network + dir_sub_quantile + 'analysis_pwm_support.np.txt')
    # # fns.append(dir_network + dir_sub_quantile + 'analysis_chip_support.bart.txt')
    # # fns.append(dir_network + dir_sub_quantile + 'analysis_pwm_support.bart.txt')
    # fns.append(dir_network + dir_sub_quantile + 'analysis_chip_support.np_bart.txt')
    # fns.append(dir_network + dir_sub_quantile + 'analysis_pwm_support.np_bart.txt')

    # # figure setup
    # colors = ['k:', 'r', 'g', 'b', 'm']
    # labels = []
    # labels.append('chance')
    # labels.append('de')
    # labels.append('lasso')
    # labels.append('np')
    # # labels.append('bart')
    # labels.append('np+bart')
    # # x_ticks = ['2k', '4k', '6k', '8k', '10k', '12k', '14k', '16k', '18k', '20k']

    # # compute chip and pwm supports
    # [eval_chip, eval_pwm] = parse_binary_gold_standard(fns, parsed.eval_method)

    # """ evaluate chip and pwm supports on binding overlap gold standard """
    # file initialization
    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_combined/'
    dir_sub_quantile = 'quantile_combine_cellCycle/analysis_binding_indep/chip.bp.np.set.sizes.top4to40k.'
    # dir_sub_quantile = 'weighed_quantile_combine_all_methods_baranski_9microarray/analysis_binding_overlap/chip.bp.np.set.sizes.top20to200k.'
    fns = []
    # fns.append(dir_network + dir_sub_quantile + 'de.txt')
    # fns.append(dir_network + dir_sub_quantile + 'lasso.txt')
    fns.append(dir_network + dir_sub_quantile + 'np.txt')
    # fns.append(dir_network + dir_sub_quantile + 'bart.txt')
    fns.append(dir_network + dir_sub_quantile + 'np_bart.txt')

    # figure setup
    colors = ['k:', 'r', 'b']
    labels = []
    labels.append('chance')
    labels.append('np')
    # labels.append('bart')
    labels.append('np+bart')
    # x_ticks = ['2k', '4k', '6k', '8k', '10k', '12k', '14k', '16k', '18k', '20k']
    x_ticks = ['4k', '8k', '12k', '16k', '20k', '24k', '28k', '32k', '36k', '40k']
    # x_ticks = ['10k', '20k', '30k', '40k', '50k', '60k', '70k', '80k', '90k', '100k']
    # x_ticks = ['20k', '40k', '60k', '80k', '100k', '120k', '140k', '160k', '180k', '200k']

    # compute chip and pwm supports
    eval_chip = [None] * (len(fns)+1)
    eval_pwm = [None] * (len(fns)+1)
    [eval_chip[0], eval_pwm[0]] = parse_chance_binding_overlap(fns[0])
    for i in range(len(fns)):
        [eval_chip[i+1], eval_pwm[i+1]] = parse_binding_overlap(fns[i], parsed.eval_method)

    """" plot figures """
    plt.figure(num=None, figsize=(15,8), dpi=80)
    plt.subplot(1,2,1)
    for i in range(len(eval_chip)):
        plt.plot(eval_chip[i], colors[i], label=labels[i])
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by ChIP')
    plt.xlim(-1, len(eval_chip[0]))
    plt.ylim(0, .5)
    # plt.legend(loc="upper right")

    plt.subplot(1,2,2)
    for i in range(len(eval_pwm)):
        plt.plot(eval_pwm[i], colors[i], label=labels[i])
    plt.xticks(range(len(eval_pwm[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by PWM')
    plt.xlim(-1, len(eval_pwm[0]))
    plt.ylim(0, .5)
    plt.legend(loc="upper right")

    plt.show()

def parse_binary_gold_standard(fns, method):
    eval_chip = numpy.zeros([len(fns)/2+1, 10])
    eval_pwm = numpy.zeros([len(fns)/2+1, 10])

    for i in range(len(fns)/2):
        chip_support = numpy.loadtxt(fns[i*2])
        pwm_support = numpy.loadtxt(fns[i*2+1]) 
        if i == 0:
            eval_chip[0,:] = chip_support[0,:]
            eval_pwm[0,:] = pwm_support[0,:]
        eval_chip[i+1,:] = chip_support[1,:]
        eval_pwm[i+1,:] = pwm_support[1,:]

    if method == 'cumulative':
        temp_eval_chip = eval_chip
        temp_eval_pwm = eval_pwm
        for i in range(10):
            eval_chip[:,i] = numpy.sum(temp_eval_chip[:,0:(i+1)], axis=1)/(i+1)
            eval_pwm[:,i] = numpy.sum(temp_eval_pwm[:,0:(i+1)], axis=1)/(i+1)
    
    return [eval_chip, eval_pwm]

def parse_binding_overlap(fn, method):
    lines = open(fn, "r").readlines()
    chip = [0] * (len(lines))
    pwm = [0] * (len(lines))

    if method == "cumulative":
        for i in range(len(lines)):
            line = lines[i].split()
            chip[i] = float(line[5])/float(line[2])
            pwm[i] = float(line[4])/float(line[2])
    
    elif method == "binned":
        for i in range(len(lines)):
            line = lines[i].split()
            if i == 0:
                chip[i] = float(line[5])/float(line[2])
                pwm[i] = float(line[4])/float(line[2])
            else:
                chip[i] = (float(line[5]) - float(prevline[5]))/(float(line[2]) - float(prevline[2]))
                pwm[i] = (float(line[4]) - float(prevline[4]))/(float(line[2]) - float(prevline[2]))
            prevline = line  

    return [chip, pwm]

def parse_chance_binding_overlap(fn):
    line = open(fn, 'r').readline()
    line = line.split()
    chip = [float(line[1])/float(line[8]) for _ in range(10)]
    pwm = [float(line[0])/float(line[7]) for _ in range(10)]

    return [chip, pwm]

if __name__ == "__main__":
    main(sys.argv)
