#!/usr/bin/python

"""
Combine mutliple tf combination of the edge scores at difference dbd cutoff.
"""

import sys
import argparse
import glob
import os.path
import numpy
import matplotlib.pyplot as plt

combination_options = ['resort', 'quantcomb']
eval_method_options = ['binned', 'cumulative']

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine chip and pwm evaluations")
    parser.add_argument('-c', '-combination', dest='combination', type=str, default='resort', help='%s' % combination_options)
    parser.add_argument('-m', '-eval_method', dest='eval_method', type=str, default='cumulative', help='%s' % eval_method_options)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    fns = []

    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/yeast_network_holstege/'
    dir_sub = 'analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.np_bart_top400k.txt')
    # fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.np_bart_top400k_tf_merged_dbd40.txt')
    # fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.np_bart_top400k_tf_merged_dbd50.txt')

    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/yeast_network_holstege_motif_incorporated/'
    dir_sub = 'fire_motifs_np_bart_bin_20/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_network_np_bart_motif_net_fire_bin_20_' + parsed.combination + '.txt')

    dir_sub = 'fire_motifs_np_bart_tf_merged_dbd50_bin_20/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_network_np_bart_motif_net_fire_bin_20_' + parsed.combination + '.txt')
    fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_network_np_bart_motif_net_fire_bin_20_tf_merged_' + parsed.combination + '.txt')
    fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_network_np_bart_tf_merged_net_fire_bin_20_tf_merged_' + parsed.combination + '.txt')

    # figure setup
    colors = ['k:', 'k', 'r', 'g', 'b', 'm']
    labels = []
    labels.append('chance')
    labels.append('np_bart')
    labels.append('np_bart + motif')
    labels.append('np_bart + tf_merged_motif')
    labels.append('np_bart + tf_merged_motif + final_net_merged')
    labels.append('np_bart_tf_merged + tf_merged_motif + final_net_merged')
    x_ticks = ['4k', '8k', '12k', '16k', '20k', '24k', '28k', '32k', '36k', '40k']

    # compute chip and pwm supports
    eval_chip = [None] * (len(fns)+1)
    eval_pwm = [None] * (len(fns)+1)
    [eval_chip[0], eval_pwm[0]] = parse_chance_binding_overlap(fns[0])
    for i in range(len(fns)):
        [eval_chip[i+1], eval_pwm[i+1]] = parse_binding_overlap(fns[i], parsed.eval_method)

    # plot figures
    plt.figure(num=None, figsize=(15,8), dpi=80)
    plt.subplot(1,2,1)
    for i in range(len(eval_chip)):
        plt.plot(eval_chip[i], colors[i], label=labels[i])
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by ChIP')
    plt.xlim(-1, len(eval_chip[0]))
    plt.ylim(0, .35)
    # plt.legend(loc="upper right")

    plt.subplot(1,2,2)
    for i in range(len(eval_pwm)):
        plt.plot(eval_pwm[i], colors[i], label=labels[i])
    plt.xticks(range(len(eval_pwm[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by PWM')
    plt.xlim(-1, len(eval_pwm[0]))
    plt.ylim(0, .35)
    plt.legend(loc="upper right")

    plt.suptitle('Yeast Holstege Network: TFs of Similar DBDs Combined for Motif Inference')

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
        temp_eval_chip = numpy.zeros([len(fns)/2+1, 10])
        temp_eval_pwm = numpy.zeros([len(fns)/2+1, 10])
        for j in range(10):
            temp_eval_chip[:,j] = numpy.sum(eval_chip[:,0:(j+1)], axis=1)/(j+1)
            temp_eval_pwm[:,j] = numpy.sum(eval_pwm[:,0:(j+1)], axis=1)/(j+1)
        eval_chip = temp_eval_chip
        eval_pwm = temp_eval_pwm

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
    chip = [float(line[1])/float(line[7]) for _ in range(10)]
    pwm = [float(line[0])/float(line[7]) for _ in range(10)]
    return [chip, pwm]

if __name__ == "__main__":
    main(sys.argv)
