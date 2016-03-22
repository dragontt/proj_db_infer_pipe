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
    parser.add_argument('-r', '-range', dest='range', default='20bins.top1.600to32k')
    parser.add_argument('-s', '-step', dest='step', default=1600, type=float)
    parser.add_argument('-c', '-combination', dest='combination', type=str, default='resort', help='%s' % combination_options)
    parser.add_argument('-m', '-eval_method', dest='eval_method', type=str, default='cumulative', help='%s' % eval_method_options)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    """ NP 1.0: different data  """
    """
    fns = []
    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/'
    # np 1.0: Hu data
    dir_sub = 'yeast_network_hu_netprophet1.0/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.combined_model.txt')
    # np 1.0: Holstege
    # dir_sub = 'yeast_network_holstege_zeke/analysis_binding_overlap/'
    # fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.zekeNp.txt')
    # # np 1.0 + bart: Holstege data
    dir_sub = 'yeast_network_holstege/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.np_bart_top400k.txt')

    # figure setup
    colors = ['k:', 'k--', 'k']
    labels = []
    labels.append('chance')
    labels.append('Netprophet 1.0: Hu data (published)')
    labels.append('Netprophet 1.0: Holstege data')
    """

    """ NP 2.0: known pwm, dbd sequences """
    
    fns = []
    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/'
    # # np 1.0: Hu data
    # dir_sub = 'yeast_network_hu_netprophet1.0/analysis_binding_overlap/'
    # fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.combined_model.txt')

    dir_sub = 'yeast_network_raw_holstege_np_global/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.combined_model_full.txt')
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.combined_model_full_tf_merged.txt')

    dir_sub = 'yeast_network_raw_holstege_bart/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.yeast_holstege_bart_full.txt')
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.yeast_holstege_bart_full_tf_merged.txt')

    dir_sub = 'yeast_network_raw_holstege_np_global/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.np_bart_combined.txt')
    dir_sub = 'yeast_network_raw_holstege_bart/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.combined_network_np_tf_merged_bart_tf_merged.txt')

    colors = ['k:', 'k', 'k--', 'g', 'g--', 'b', 'b--']
    labels = ["chance", "np", "np_tf_merged", "bart", "bart_tf_merged", "np + bart", "np_tf_merged + bart_tf_merged"]
    

    """ Figure setup """
    
    dir_figures = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/yeast_analysis_results/'
    # figure_title = 'Yeast Holstege Network: TF Merging and Combination with Motifs'

    # compute chip and pwm supports
    eval_chip = [None] * (len(fns)+1)
    eval_pwm = [None] * (len(fns)+1)
    [eval_chip[0], eval_pwm[0]] = parse_chance_binding_overlap(fns[1])
    for i in range(len(fns)):
        [eval_chip[i+1], eval_pwm[i+1]] = parse_binding_overlap(fns[i], parsed.eval_method)
    print 'chip chance:', eval_chip[0][0], 'pwm chance:', eval_pwm[0][0]

    x_ticks = [format(float(i)*parsed.step/320, '.0f') for i in range(1,len(eval_chip[0])+1)]
    

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
    plt.xlabel('Average number of predicted targets per TF in the genome')
    plt.ylabel('Interactions supported by ChIP (%)')
    plt.xlim(-1, len(eval_chip[0])+1)
    plt.ylim(2, 50)
    plt.yticks(numpy.arange(0,51,5))
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)

    ax = plt.subplot(1,2,2)
    for i in range(len(eval_pwm)):
        ax.plot(eval_pwm[i], colors[i], label=labels[i], linewidth=2.0)
    ax.scatter(18, 10, s=75, c='c')
    ax.annotate('ChIP network', xy=(17.75, 9.75), xycoords='data', xytext=(-50, -30), textcoords='offset points', arrowprops=dict(arrowstyle="->"))
    plt.xticks(range(len(eval_pwm[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF in the genome')
    plt.ylabel('Interactions supported by PWM (%)')
    plt.xlim(-1, len(eval_pwm[0])+1)
    plt.ylim(5, 30)
    plt.yticks(numpy.arange(0,31,5))
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., scatterpoints=1)

    plt.savefig(dir_figures + parsed.figure_name + '.pdf', fmt='pdf')
    

    """ Plot ChIP supoort """
    """
    # broken axis
    ylim1, ylim2 = [10,50], [2.5,4.5]
    ylim1_ratio = (ylim1[1]-ylim1[0])/(ylim2[1]-ylim2[0]+ylim1[1]-ylim1[0])
    ylim2_ratio = (ylim2[1]-ylim2[0])/(ylim2[1]-ylim2[0]+ylim1[1]-ylim1[0])

    gs = gridspec.GridSpec(2,1,height_ratios=[ylim1_ratio, ylim2_ratio])
    fig = plt.figure(num=None, figsize=(6,6), dpi=150)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    for i in range(len(eval_chip)):
        ax1.plot(eval_chip[i], colors[i], label=labels[i], linewidth=2.0)
        ax2.plot(eval_chip[i], colors[i], label=labels[i], linewidth=2.0)
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF in the genome')
    plt.ylabel('Interactions supported by ChIP (%)')
    plt.xlim(-1, len(eval_chip[0]))
    ax2.yaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)
    # plt.ylim(13.5,37)
    plt.subplots_adjust(hspace=.1)
    # ax1.legend(loc="upper right")
    
    # hide the spines between ax and ax2
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop='off') # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(color='k', clip_on=False)
    xlim = ax2.get_xlim()
    dx = .02*(xlim[1]-xlim[0])
    dy = .01*(ylim1[1]-ylim1[0])/ylim1_ratio
    ax1.plot((xlim[0]-dx,xlim[0]+dx), (ylim1[0]-dy,ylim1[0]+dy), **kwargs)
    ax1.plot((xlim[1]-dx,xlim[1]+dx), (ylim1[0]-dy,ylim1[0]+dy), **kwargs)
    dy = .01*(ylim2[1]-ylim2[0])/ylim2_ratio
    ax2.plot((xlim[0]-dx,xlim[0]+dx), (ylim2[1]-dy,ylim2[1]+dy), **kwargs)
    ax2.plot((xlim[1]-dx,xlim[1]+dx), (ylim2[1]-dy,ylim2[1]+dy), **kwargs)

    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)
    ax1.set_ylim(ylim1)
    ax2.set_ylim(ylim2)
    for label in ax2.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    ax2.get_yaxis().set_ticks([])

    # ax.set_xscale('log')
    # plt.xlim(0, len(eval_chip[0]))

    plt.savefig(dir_figures + parsed.figure_name + '_chip.pdf', fmt='pdf')
    """

    """ Plot PWM supoort """
    """
    # broken axis
    ylim1, ylim2 = [12,28], [5.5,7]
    ylim1_ratio = (ylim1[1]-ylim1[0])/(ylim2[1]-ylim2[0]+ylim1[1]-ylim1[0])
    ylim2_ratio = (ylim2[1]-ylim2[0])/(ylim2[1]-ylim2[0]+ylim1[1]-ylim1[0])

    gs = gridspec.GridSpec(2,1,height_ratios=[ylim1_ratio, ylim2_ratio])
    fig = plt.figure(num=None, figsize=(6,6), dpi=150)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    for i in range(len(eval_chip)):
        ax1.plot(eval_pwm[i], colors[i], label=labels[i], linewidth=2.0)
        ax2.plot(eval_pwm[i], colors[i], label=labels[i], linewidth=2.0)
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF in the genome')
    plt.ylabel('Interactions supported by PWM (%)')
    plt.xlim(-1, len(eval_chip[0]))
    ax2.yaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)
    # plt.ylim(13.5,37)
    plt.subplots_adjust(hspace=.1)
    ax1.legend(loc="upper right")
    
    # hide the spines between ax and ax2
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop='off') # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(color='k', clip_on=False)
    xlim = ax2.get_xlim()
    dx = .02*(xlim[1]-xlim[0])
    dy = .01*(ylim1[1]-ylim1[0])/ylim1_ratio
    ax1.plot((xlim[0]-dx,xlim[0]+dx), (ylim1[0]-dy,ylim1[0]+dy), **kwargs)
    ax1.plot((xlim[1]-dx,xlim[1]+dx), (ylim1[0]-dy,ylim1[0]+dy), **kwargs)
    dy = .01*(ylim2[1]-ylim2[0])/ylim2_ratio
    ax2.plot((xlim[0]-dx,xlim[0]+dx), (ylim2[1]-dy,ylim2[1]+dy), **kwargs)
    ax2.plot((xlim[1]-dx,xlim[1]+dx), (ylim2[1]-dy,ylim2[1]+dy), **kwargs)

    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)
    ax1.set_ylim(ylim1)
    ax2.set_ylim(ylim2)
    for label in ax2.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    ax2.get_yaxis().set_ticks([])

    # ax.set_xscale('log')
    # plt.xlim(0, len(eval_pwm[0]))

    plt.savefig(dir_figures + parsed.figure_name + '_pwm.pdf', fmt='pdf')
    """

    """ Bar plot at specific targets per tf level """
    """
    fns = []
    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/'
    # np 1.0 + bart: Holstege data
    dir_sub = 'yeast_network_holstege/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.np_bart_top400k.txt')
    # np + tf_merging
    dir_sub = 'yeast_network_holstege/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.np_bart_top400k_tf_merged_dbd50.txt')
    # np + cis-bp motifs
    dir_sub = 'yeast_network_holstege_motif_incorporated/cisbp_dbd_only_holstege_np_bart_tf_merged_dbd50/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.combined_network_np_motif_net_dbd_only_tf_merged_dbd50.txt')
    # np + fire
    dir_sub = 'yeast_network_holstege_motif_incorporated/fire_motifs_np_bart_tf_merged_dbd50_bin_20/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.combined_network_np_bart_tf_merged_net_fire_bin_20_tf_merged_resort.txt')
    # np + fire_ortho
    dir_sub = 'yeast_network_holstege_motif_incorporated/fire_ortho_scer+spar_motifs_bin_20/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.combined_np_bart_tf_merged_motif_net_tf_merged.txt')
     # np + known pwm
    dir_sub = 'yeast_network_holstege_motif_incorporated/scertf_known_motif/analysis_binding_overlap/'
    fns.append(dir_network + dir_sub + 'analysis.'+ parsed.range +'.combined_np_bart_tf_merged_motif_net_tf_merged.txt')

    # figure setup
    colors = ['0.5', 'k', 'r', 'b', 'm', 'c', 'g']
    labels = []
    labels.append('chance')
    labels.append('NP 1.0: Holstege data')
    labels.append('NP 1.0 + weighted averaging (WA)')
    labels.append('NP 1.0 + WA + CIS-BP motifs + WA')
    labels.append('NP 1.0 + WA + FIRE motifs\n(target species only) + WA')
    labels.append('NP 1.0 + WA + FIRE motifs\n(target + orthologs) + WA')
    labels.append('NP 1.0 + WA + known motifs + WA')

    # compute chip and pwm supports
    eval_chip = [None] * (len(fns)+1)
    eval_pwm = [None] * (len(fns)+1)
    [eval_chip[0], eval_pwm[0]] = parse_chance_binding_overlap(fns[0])
    for i in range(len(fns)):
        [eval_chip[i+1], eval_pwm[i+1]] = parse_binding_overlap(fns[i], parsed.eval_method)
    eval_chip = numpy.array(eval_chip[::-1])
    eval_pwm = numpy.array(eval_pwm[::-1])
    colors = colors[::-1]
    labels = labels[::-1]

    # plot figure
    dir_figures = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/yeast_analysis_results/'

    xlevel = 25
    xlevel_index = xlevel/5-1 

    fig = plt.figure(num=None, figsize=(20,5), dpi=150)
    fig.subplots_adjust(wspace=.5)

    y_pos = numpy.arange(len(labels))

    plt.subplot(1,2,1)
    for i in y_pos:
        plt.barh(y_pos[i], eval_chip[:,xlevel_index][i], align='center', alpha=1, color=colors[i], edgecolor = "white")
    plt.ylabel('Network mapping procedures')
    plt.xlabel('Interactions supported by ChIP (%)')
    plt.xlim([0,28])
    plt.xticks(numpy.arange(0,29,2))
    plt.tick_params(labelleft='off')
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()

    ax = plt.subplot(1,2,2)
    for i in y_pos:
        plt.barh(y_pos[i], eval_pwm[:,xlevel_index][i], align='center', alpha=1, color=colors[i], edgecolor = "white")
    plt.yticks(y_pos, labels)
    ax.set_yticklabels(labels, horizontalalignment='center', position=(-.25,.5))
    plt.xlabel('Interactions supported by PWM (%)')
    plt.xlim([0,22])
    plt.xticks(numpy.arange(0,23,2))
    plt.gca().invert_yaxis()

    # plt.savefig(dir_figures + parsed.figure_name + '_bar_' + str(xlevel) +'.pdf' , fmt='pdf')
    plt.savefig(dir_figures + parsed.figure_name +'.pdf' , fmt='pdf')
    """

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
