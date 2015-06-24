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
    parser.add_argument('-e', '-eval_method', dest='eval_method', type=str, default='cumulative')
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    """ evaluate chip and pwm supports on binary gold standard """
    # file initialization
    fns = []

    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_global_shrinkage/'
    dir_sub = 'analysis_flynet_top4to40k/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.combined_model_full.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.combined_model_full.txt')

    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_motif_incorporated/'
    dir_sub = 'cisbp_-2000_+200_fimo_dbd_cutoff_cellCycle_np/analysis_flynet_top4to40k/'

    fns.append(dir_network + dir_sub + 'analysis_chip_support.combined_network_np_motif_net_dbd_cutoff_40.0_geomean.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.combined_network_np_motif_net_dbd_cutoff_40.0_geomean.txt')
    fns.append(dir_network + dir_sub + 'analysis_chip_support.combined_network_np_motif_net_dbd_cutoff_50.0_geomean.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.combined_network_np_motif_net_dbd_cutoff_50.0_geomean.txt')
    fns.append(dir_network + dir_sub + 'analysis_chip_support.combined_network_np_motif_net_dbd_cutoff_60.0_geomean.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.combined_network_np_motif_net_dbd_cutoff_60.0_geomean.txt')

    fns.append(dir_network + dir_sub + 'analysis_chip_support.combined_network_np_motif_net_dbd_cutoff_40.0_resort.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.combined_network_np_motif_net_dbd_cutoff_40.0_resort.txt')
    fns.append(dir_network + dir_sub + 'analysis_chip_support.combined_network_np_motif_net_dbd_cutoff_50.0_resort.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.combined_network_np_motif_net_dbd_cutoff_50.0_resort.txt')
    fns.append(dir_network + dir_sub + 'analysis_chip_support.combined_network_np_motif_net_dbd_cutoff_60.0_resort.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.combined_network_np_motif_net_dbd_cutoff_60.0_resort.txt')

    fns.append(dir_network + dir_sub + 'analysis_chip_support.combined_network_np_motif_net_dbd_cutoff_40.0_quantcomb.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.combined_network_np_motif_net_dbd_cutoff_40.0_quantcomb.txt')
    fns.append(dir_network + dir_sub + 'analysis_chip_support.combined_network_np_motif_net_dbd_cutoff_50.0_quantcomb.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.combined_network_np_motif_net_dbd_cutoff_50.0_quantcomb.txt')
    fns.append(dir_network + dir_sub + 'analysis_chip_support.combined_network_np_motif_net_dbd_cutoff_60.0_quantcomb.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.combined_network_np_motif_net_dbd_cutoff_60.0_quantcomb.txt')

    # figure setup
    colors = ['k:', 'k', 'r:', 'g:', 'b:', 'r--', 'g--', 'b--', 'r', 'g', 'b']
    labels = []
    labels.append('chance')
    labels.append('np')
    labels.append('np + motif_cutoff_40_geomean')
    labels.append('np + motif_cutoff_50_geomean')
    labels.append('np + motif_cutoff_60_geomean')
    labels.append('np + motif_cutoff_40_resort')
    labels.append('np + motif_cutoff_50_resort')
    labels.append('np + motif_cutoff_60_resort')
    labels.append('np + motif_cutoff_40_quantcomb')
    labels.append('np + motif_cutoff_50_quantcomb')
    labels.append('np + motif_cutoff_60_quantcomb')
    x_ticks = ['4k', '8k', '12k', '16k', '20k', '24k', '28k', '32k', '36k', '40k']

    # compute chip and pwm supports
    [eval_chip, eval_pwm] = parse_binary_gold_standard(fns, parsed.eval_method)

    """ evaluate chip and pwm supports on binding overlap gold standard """
    # # file initialization
    # fns = []

    # dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_global_shrinkage/'
    # dir_sub = 'analysis_binding_overlap/'
    # fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_model_full.txt')

    # dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_motif_incorporated/'
    # dir_sub = 'cisbp_-2000_+200_fimo_dbd_cutoff_cellCycle_np/analysis_binding_overlap/'
    # fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_network_np_motif_net_dbd_cutoff_40.0.txt')
    # fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_network_np_motif_net_dbd_cutoff_50.0.txt')
    # fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_network_np_motif_net_dbd_cutoff_60.0.txt')

    # dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_motif_incorporated/'
    # dir_sub = 'cisbp_-2000_+200_fimo_dbd_knn_cellCycle_np/analysis_binding_overlap/'
    # fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_network_np_motif_net_dbd_knn_10.txt')
    # fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_network_np_motif_net_dbd_knn_20.txt')
    # fns.append(dir_network + dir_sub + 'chip.bp.np.set.sizes.top4to40k.combined_network_np_motif_net_dbd_knn_30.txt')

    # # figure setup
    # colors = ['k:', 'k', 'r', 'g', 'b', 'r--', 'g--', 'b--']
    # labels = []
    # labels.append('chance')
    # labels.append('np')
    # labels.append('np + motif_cutoff_40')
    # labels.append('np + motif_cutoff_50')
    # labels.append('np + motif_cutoff_60')
    # labels.append('np + motif_knn_10')
    # labels.append('np + motif_knn_20')
    # labels.append('np + motif_knn_30')
    # x_ticks = ['4k', '8k', '12k', '16k', '20k', '24k', '28k', '32k', '36k', '40k']
    # # x_ticks = ['10k', '20k', '30k', '40k', '50k', '60k', '70k', '80k', '90k', '100k']
    # # x_ticks = ['20k', '40k', '60k', '80k', '100k', '120k', '140k', '160k', '180k', '200k']

    # # compute chip and pwm supports
    # eval_chip = [None] * (len(fns)+1)
    # eval_pwm = [None] * (len(fns)+1)
    # [eval_chip[0], eval_pwm[0]] = parse_chance_binding_overlap(fns[0])
    # for i in range(len(fns)):
    #     [eval_chip[i+1], eval_pwm[i+1]] = parse_binding_overlap(fns[i], parsed.eval_method)

    """" plot figures """
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
    chip = [float(line[1])/float(line[8]) for _ in range(10)]
    pwm = [float(line[0])/float(line[7]) for _ in range(10)]

    return [chip, pwm]

if __name__ == "__main__":
    main(sys.argv)
