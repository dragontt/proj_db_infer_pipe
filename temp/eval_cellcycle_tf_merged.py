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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine chip and pwm evaluations")
    parser.add_argument('-n', '-figure_name', dest='figure_name', default='temp_result')
    parser.add_argument('-r', '-range', dest='range', default='20bins.top4.845to96.9k')
    parser.add_argument('-s', '-step', dest='step', default=4845, type=float)
    parser.add_argument('-c', '-combination', dest='combination', type=str, default='resort')
    parser.add_argument('-m', '-eval_method', dest='eval_method', type=str, default='cumulative')
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)
    
    fns = []
    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/'

    # np
    dir_sub = 'fly_network_cellCycle_np_bart/analysis_compiled_chip_flynet_pwm/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.np_full.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.np_full.txt')

    # np_tf_merged + bart_tf_merged
    dir_sub = 'fly_network_cellCycle_np_bart/analysis_compiled_chip_flynet_pwm/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_np_tf_merged_bart_tf_merged.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_np_tf_merged_bart_tf_merged.txt')

    # np + bart + known motifs
    dir_sub = 'fly_network_cellCycle_np_bart_motif_incorporated/known_cisbp/analysis_compiled_chip_flynet_pwm/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_np_bart_motif_tf_merged.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_np_bart_motif_tf_merged.txt')

    # figure setup
    colors = ['k:', 'k', 'b', 'g']
    labels = []
    labels.append('chance')
    labels.append('NP')
    labels.append('NP + BART')
    labels.append('NP + BART + known_motifs')

    """ Figure setup """
    
    dir_figures = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_analysis_results/'
    # figure_title = 'Fly CellCycle Network: Combination with Inferred Motifs'

    # compute chip and pwm supports
    [eval_chip, eval_pwm] = parse_binary_gold_standard(fns, parsed.eval_method)
    print 'chip chance:', eval_chip[0][0], 'pwm chance:', eval_pwm[0][0]

    # targets per tf
    x_ticks = [format(float(i)*parsed.step/969, '.0f') for i in range(1,len(eval_chip[0])+1)]

    """ Regular support plot """
    
    # plot figures
    # fig = plt.figure(num=None, figsize=(25,8), dpi=80)
    # fig.subplots_adjust(right=.7)

    fig = plt.figure(num=None, figsize=(18,6), dpi=80)
    fig.subplots_adjust(right=.7)

    ax = plt.subplot(1,2,1)
    for i in range(len(eval_chip)):
        ax.plot(eval_chip[i], colors[i], label=labels[i], linewidth=1.5)
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF in the genome')
    plt.ylabel('Interactions supported by ChIP (%)')
    plt.xlim(-1, len(eval_chip[0])+1)
    # plt.ylim(0,60)
    # plt.yticks(numpy.arange(0,60.5,10))
    plt.ylim(0,90)
    plt.yticks(numpy.arange(0,90.5,10))
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)

    ax = plt.subplot(1,2,2)
    for i in range(len(eval_pwm)):
        ax.plot(eval_pwm[i], colors[i], label=labels[i], linewidth=1.5)
    # ax.scatter(11,12.4, s=75, c='c')
    # ax.annotate('ChIP network', xy=(10.75, 12), xycoords='data', xytext=(-50, -30), textcoords='offset points', arrowprops=dict(arrowstyle="->"))
    plt.xticks(range(len(eval_pwm[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF in the genome')
    plt.ylabel('Interactions supported by PWM (%)')
    plt.xlim(-1, len(eval_pwm[0])+1)
    plt.ylim(0,12)
    plt.yticks(numpy.arange(0,12.5,2))
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., scatterpoints=1)

    plt.savefig(dir_figures + parsed.figure_name + '.pdf', fmt='pdf')


def parse_binary_gold_standard(fns, method):
    evalPoints = numpy.loadtxt(fns[0]).shape[1]

    eval_chip = numpy.zeros([len(fns)/2+1, evalPoints])
    eval_pwm = numpy.zeros([len(fns)/2+1, evalPoints])

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
            for j in range(1,evalPoints):
                eval_chip[i+1,j] = (chip_support[3,j]-chip_support[3,j-1])/(chip_support[2,j]-chip_support[2,j-1])
                eval_pwm[i+1,j] = (pwm_support[3,j]-pwm_support[3,j-1])/(pwm_support[2,j]-pwm_support[2,j-1])

    return [eval_chip*100, eval_pwm*100]

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
