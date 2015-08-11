#!/usr/bin/python

"""
Combine mutliple tf combination of the edge scores at difference dbd cutoff.

eval_type: flynet, binding_indep
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
    parser.add_argument('-t', '-eval_type', dest='eval_type', type=str, default=
        'flynet')
    parser.add_argument('-c', '-combination', dest='combination', type=str, default='resort')
    parser.add_argument('-m', '-eval_method', dest='eval_method', type=str, default='cumulative')
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    if parsed.eval_type == "flynet":

        # evaluate chip and pwm supports on binary gold standard 
        # file initialization
        fns = []

        dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_combined_baranski_cellcycle/'
        dir_sub = 'analysis_compiled_chip_flynet_pwm/'
        fns.append(dir_network + dir_sub + 'analysis_chip_support.top12to120k.np_combined.txt')
        fns.append(dir_network + dir_sub + 'analysis_pwm_support.top12to120k.np_combined.txt')
        fns.append(dir_network + dir_sub + 'analysis_chip_support.top12to120k.np_combined_tf_merged_dbd40.txt')
        fns.append(dir_network + dir_sub + 'analysis_pwm_support.top12to120k.np_combined_tf_merged_dbd40.txt')
        fns.append(dir_network + dir_sub + 'analysis_chip_support.top12to120k.np_combined_tf_merged_dbd50.txt')
        fns.append(dir_network + dir_sub + 'analysis_pwm_support.top12to120k.np_combined_tf_merged_dbd50.txt')

        # figure setup
        colors = ['k:', 'k', 'r', 'g']
        labels = []
        labels.append('chance')
        labels.append('np_combined')
        labels.append('np_combined_tf_merged_dbd40')
        labels.append('np_combined_tf_merged_dbd50')

        # x_ticks = ['4k', '8k', '12k', '16k', '20k', '24k', '28k', '32k', '36k', '40k']
        x_ticks = ['12k', '24k', '36k', '48k', '60k', '72k', '84k', '96k', '108k', '120k']

        # compute chip and pwm supports
        [eval_chip, eval_pwm] = parse_binary_gold_standard(fns, parsed.eval_method)

    elif parsed.eval_type == "binding_indep":
        # evaluate chip and pwm supports on binding overlap gold standard
        # file initialization
        fns = []

        dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_combined_baranski_cellcycle/'
        dir_sub = 'analysis_binding_indep/'
        fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_model_full.txt')

        # figure setup
        colors = ['k:', 'k', 'r', 'g', 'g--', 'b', 'b--']
        labels = []
        labels.append('chance')
        labels.append('np')
        labels.append('np + known_motif')
        labels.append('np + fire_motif')
        labels.append('np + fire_motif_tf_merged')
        labels.append('np + database_motif')
        labels.append('np + database_motif_tf_merged')
        x_ticks = ['4k', '8k', '12k', '16k', '20k', '24k', '28k', '32k', '36k', '40k']

        # compute chip and pwm supports
        eval_chip = [None] * (len(fns)+1)
        eval_pwm = [None] * (len(fns)+1)
        [eval_chip[0], eval_pwm[0]] = parse_chance_binding_indep(fns[0])
        for i in range(len(fns)):
            [eval_chip[i+1], eval_pwm[i+1]] = parse_binding_info(fns[i], parsed.eval_method)

    # plot figures
    plt.figure(num=None, figsize=(15,8), dpi=80)
    plt.subplot(1,2,1)
    for i in range(len(eval_chip)):
        plt.plot(eval_chip[i], colors[i], label=labels[i])
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by ChIP')
    plt.xlim(-1, len(eval_chip[0]))
    plt.ylim(0, .6)
    # plt.legend(loc="upper right")

    plt.subplot(1,2,2)
    for i in range(len(eval_pwm)):
        plt.plot(eval_pwm[i], colors[i], label=labels[i])
    plt.xticks(range(len(eval_pwm[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by PWM')
    plt.xlim(-1, len(eval_pwm[0]))
    plt.ylim(0, .3)
    plt.legend(loc="upper right")

    # plt.suptitle(parsed.eval_type + '_tf_merged_dbd50_' + parsed.combination)

    plt.show()

def parse_binary_gold_standard(fns, method):
    eval_chip = numpy.zeros([len(fns)/2+1, 10])
    eval_pwm = numpy.zeros([len(fns)/2+1, 10])

    for i in range(len(fns)/2):
        chip_support = numpy.loadtxt(fns[i*2], dtype=str)
        temp_index = numpy.where(chip_support == 'NA')
        chip_support[temp_index] = '0'
        chip_support = numpy.array(chip_support, dtype=float)

        pwm_support = numpy.loadtxt(fns[i*2+1], dtype=str) 
        temp_index = numpy.where(pwm_support == 'NA')
        pwm_support[temp_index] = '0'
        pwm_support = numpy.array(pwm_support, dtype=float)
        
        if i == 0:
            eval_chip[0,:] = chip_support[1,:]/chip_support[0,:]
            eval_pwm[0,:] = pwm_support[1,:]/pwm_support[0,:]

        if method == 'cumulative':
            eval_chip[i+1,:] = chip_support[3,:]/chip_support[2,:]
            eval_pwm[i+1,:] = pwm_support[3,:]/pwm_support[2,:]

        elif method == 'binned':
            eval_chip[i+1,0] = chip_support[3,0]/chip_support[2,0]
            eval_pwm[i+1,0] = pwm_support[3,0]/pwm_support[2,0]
            for j in range(1,10):
                eval_chip[i+1,j] = (chip_support[3,j]-chip_support[3,j-1])/(chip_support[2,j]-chip_support[2,j-1])
                eval_pwm[i+1,j] = (pwm_support[3,j]-pwm_support[3,j-1])/(pwm_support[2,j]-pwm_support[2,j-1])

    return [eval_chip, eval_pwm]

def parse_binding_info(fn, method):
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

def parse_chance_binding_indep(fn):
    line = open(fn, 'r').readline()
    line = line.split()
    chip = [float(line[1])/float(line[8]) for _ in range(10)]
    pwm = [float(line[0])/float(line[7]) for _ in range(10)]
    return [chip, pwm]

if __name__ == "__main__":
    main(sys.argv)
