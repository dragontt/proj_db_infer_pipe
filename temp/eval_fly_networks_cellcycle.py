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
import matplotlib.gridspec as gridspec

combination_options = ['resort', 'quantcomb']
eval_method_options = ['binned', 'cumulative']

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine chip and pwm evaluations")
    parser.add_argument('-n', '-figure_name', dest='figure_name', default='temp_result')
    parser.add_argument('-r', '-range', dest='range', default='top4to40k')
    parser.add_argument('-s', '-step', dest='step', default=1000, type=float)
    parser.add_argument('-m', '-eval_method', dest='eval_method', type=str, default='cumulative', help='%s' % eval_method_options)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    """ Parse evaluation data """

    fns = []

    # original network
    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_analysis_results/networks_cellcycle/analysis_binding_overlap/'
    # dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_analysis_results/networks_cellcycle/analysis_binding_overlap_low_stringency/'
    fns.append(dir_network + 'analysis.'+ parsed.range +'.np.txt')
    fns.append(dir_network + 'analysis.'+ parsed.range +'.np_tf_merged.txt')
    fns.append(dir_network + 'analysis.'+ parsed.range +'.np_known_motif.txt')
    fns.append(dir_network + 'analysis.'+ parsed.range +'.np_cisbp_dbd_only.txt')
    fns.append(dir_network + 'analysis.'+ parsed.range +'.np_fire_dmel.txt')
    fns.append(dir_network + 'analysis.'+ parsed.range +'.np_fire_dmel_dsim_dsec.txt')

    """ Figure setup """
    
    dir_figures = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_analysis_results/'
    figure_title = 'Yeast Holstege Network: TF Merging and Combination with Motifs'

    # figure setup
    colors = ['k:', 'k', 'c', 'r', 'g', 'b', 'm']
    labels = []
    labels.append('chance')
    labels.append('NetProphet (NP) 1.0')
    labels.append('NP 1.0 + Weighted averaging (WA)')
    labels.append('NP 1.0 + WA + Known motifs')
    labels.append('NP 1.0 + WA + Database-inferred motifs')
    labels.append('NP 1.0 + WA + FIRE-inferred motifs')
    labels.append('NP 1.0 + WA + FIRE-inferred motifs with orthologs')

    # compute chip and pwm supports
    eval_chip = [None] * (len(fns)+1)
    eval_pwm = [None] * (len(fns)+1)
    [eval_chip[0], eval_pwm[0]] = parse_chance_binding_overlap(fns[0])
    for i in range(len(fns)):
        [eval_chip[i+1], eval_pwm[i+1]] = parse_binding_overlap(fns[i], parsed.eval_method)
    print 'chip chance:', eval_chip[0][0], 'pwm chance:', eval_pwm[0][0]

    x_ticks = [format(float(i)*parsed.step/969, '.0f') for i in range(1,len(eval_chip[0])+1)]


    """ Regular support plot """
    # plot figures
    # fig = plt.figure(num=None, figsize=(25,8), dpi=80)
    # fig.subplots_adjust(right=.7)

    fig = plt.figure(num=None, figsize=(18,6), dpi=80)
    fig.subplots_adjust(right=.7)

    ax = plt.subplot(1,2,1)
    for i in range(len(eval_chip)):
        ax.plot(eval_chip[i], colors[i], label=labels[i], linewidth=2.0)
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF genome')
    plt.ylabel('Interactions supported by ChIP (%)')
    plt.xlim(-1, len(eval_chip[0]))
    plt.ylim(2, 50)
    plt.yticks(numpy.arange(0,51,5))
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)

    ax = plt.subplot(1,2,2)
    for i in range(len(eval_pwm)):
        ax.plot(eval_pwm[i], colors[i], label=labels[i], linewidth=2.0)
    ax.scatter(18, 10, s=75, c='c', label='ChIP network')
    plt.xticks(range(len(eval_pwm[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF genome')
    plt.ylabel('Interactions supported by PWM (%)')
    plt.xlim(-1, len(eval_pwm[0]))
    plt.ylim(5, 30)
    plt.yticks(numpy.arange(0,31,5))
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    # plt.legend(loc="upper right", scatterpoints=1)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., scatterpoints=1)

    plt.savefig(dir_figures + parsed.figure_name + '.pdf', fmt='pdf')

    
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

    return [eval_chip*100, eval_pwm*100]

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
    return [numpy.array(chip)*100, numpy.array(pwm)*100]

def parse_chance_binding_overlap(fn):
    evalPoints = numpy.loadtxt(fn).shape[0]

    line = open(fn, 'r').readline()
    line = line.split()
    chip = [float(line[1])/float(line[7]) for _ in range(evalPoints)]
    pwm = [float(line[0])/float(line[7]) for _ in range(evalPoints)]
    return [numpy.array(chip)*100, numpy.array(pwm)*100]

if __name__ == "__main__":
    main(sys.argv)
