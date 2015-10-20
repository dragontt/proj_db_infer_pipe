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
import matplotlib.gridspec as gridspec

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine chip and pwm evaluations")
    parser.add_argument('-n', '-figure_name', dest='figure_name', default='temp_result')
    parser.add_argument('-r', '-range', dest='range', default='top4to40k')
    parser.add_argument('-s', '-step', dest='step', default=1000, type=float)
    parser.add_argument('-c', '-combination', dest='combination', type=str, default='resort')
    parser.add_argument('-m', '-eval_method', dest='eval_method', type=str, default='cumulative')
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    """ Parse evaluation data """

    # evaluate chip and pwm supports on binary gold standard 
    # file initialization
    fns = []

    # np original
    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_global_shrinkage/'
    dir_sub = 'analysis_compiled_chip_flynet_pwm/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_model_full.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_model_full.txt')
    # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_model_full_tf_merged_pid50.txt')
    # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_model_full_tf_merged_pid50.txt')

    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_motif_incorporated/'

    # fire inference
    dir_sub = 'fire_motifs_np_bin_20/analysis_compiled_chip_flynet_pwm/'
    # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_motif_net_fire_np_bin_20_resort.txt')
    # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_motif_net_fire_np_bin_20_resort.txt')

    dir_sub = 'fire_motifs_np_tf_merged_dbd50_bin_20/analysis_compiled_chip_flynet_pwm/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_motif_net_fire_np_bin_20_resort_tf_merged.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_motif_net_fire_np_bin_20_resort_tf_merged.txt')

    dir_sub = 'fire_ortho_motifs_np_bin_20/analysis_compiled_chip_flynet_pwm/'
    # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')
    # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')
    dir_sub = 'fire_ortho_dmel+Dsim+Dsec_motifs_np_bin_20/analysis_compiled_chip_flynet_pwm/'
    # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')
    # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')
    # dir_sub = 'fire_ortho_dmel+Dsim+Dsec+Dyak+Dere+Dana_motifs_np_bin_20/analysis_compiled_chip_flynet_pwm/'
    # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')
    # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')

    # cisbp inference
    dir_sub = 'cisbp_-2000_+200_fimo_dbd_cutoff_cellCycle_np/analysis_compiled_chip_flynet_pwm/'
    # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_motif_net_dbd_cutoff_40.0_resort_0.txt')
    # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_motif_net_dbd_cutoff_40.0_resort_0.txt')

    dir_sub = 'cisbp_-2000_+200_fimo_dbd_cutoff_cellCycle_np_tf_merged_dbd50/analysis_compiled_chip_flynet_pwm/'
    # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_motif_net_dbd_cutoff_40_resort_tf_merged.txt')
    # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_motif_net_dbd_cutoff_40_resort_tf_merged.txt')

    # known motifs
    dir_sub = 'cisbp_-2000_+200_fimo_known_motif/analysis_compiled_chip_flynet_pwm/'
    # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_motif_net_known_motif_resort.txt')
    # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_motif_net_known_motif_resort.txt')
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_tf_merged_motif_net_known_motif_tf_merged_resort.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_tf_merged_motif_net_known_motif_tf_merged_resort.txt')


    """ Figure setup """

    dir_figures = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_analysis_results/'
    figure_title = 'Fly CellCycle Network: Combination with Inferred Motifs'

    # figure setup
    # colors = ['k:', 'k', 'k--', 'r', 'r--', 'g', 'g--', 'g:', 'b', 'b--', 'm', 'm--' ]
    colors = ['k:', 'k--', 'r', 'g--']
    # colors = ['k:', 'k--', 'k', 'r', 'g', 'b', 'm--']
    labels = []
    labels.append('chance')
    labels.append('np')
    # labels.append('np + tf_score_summing')
    # labels.append('np + fire_dmel_motifs')
    # labels.append('np + tf_score_summing + fire_motifs_dmel')
    labels.append('np + tf_score_summing + inferred_motifs')
    # labels.append('np + tf_score_summing + fire_indiv_ortho_motifs')
    # labels.append('np + tf_score_summing + fire_motifs_ortho_dmel+dsim+dsec')
    # labels.append('np + tf_score_summing + fire_comb_ortho_motifs _5species')
    # labels.append('np + cisbp_motifs')
    # labels.append('np + tf_score_summing + cisbp_motifs')
    # labels.append('np + known_motifs')
    labels.append('np + tf_score_summing + known_motifs')

    # compute chip and pwm supports
    [eval_chip, eval_pwm] = parse_binary_gold_standard(fns, parsed.eval_method)
    print 'chip chance:', eval_chip[0][0], 'pwm chance:', eval_pwm[0][0]

    # targets per tf
    x_ticks = [format(float(i)*parsed.step/969, '.0f') for i in range(1,len(eval_chip[0]+1))]

    """ Plot ChIP supoort """
    # broken axis
    ylim1, ylim2 = [25,52], [14.,16.]
    ylim1_ratio = (ylim1[1]-ylim1[0])/(ylim2[1]-ylim2[0]+ylim1[1]-ylim1[0])
    ylim2_ratio = (ylim2[1]-ylim2[0])/(ylim2[1]-ylim2[0]+ylim1[1]-ylim1[0])

    gs = gridspec.GridSpec(2,1,height_ratios=[ylim1_ratio, ylim2_ratio])
    fig = plt.figure(num=None, figsize=(8,8), dpi=80)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    for i in range(len(eval_chip)):
        ax1.plot(eval_chip[i], colors[i], label=labels[i], linewidth=2.0)
        ax2.plot(eval_chip[i], colors[i], label=labels[i], linewidth=2.0)
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF genome')
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
    ax2.get_yaxis().set_ticks([])

    # ax.set_xscale('log')
    # plt.xlim(0, len(eval_chip[0]))

    plt.savefig(dir_figures + parsed.figure_name + '_chip.pdf', fmt='pdf')

    """ Plot PWM supoort """
    # broken axis
    ylim1, ylim2 = [6,13], [4.85,5.25]
    ylim1_ratio = (ylim1[1]-ylim1[0])/(ylim2[1]-ylim2[0]+ylim1[1]-ylim1[0])
    ylim2_ratio = (ylim2[1]-ylim2[0])/(ylim2[1]-ylim2[0]+ylim1[1]-ylim1[0])

    gs = gridspec.GridSpec(2,1,height_ratios=[ylim1_ratio, ylim2_ratio])
    fig = plt.figure(num=None, figsize=(8,8), dpi=80)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    for i in range(len(eval_chip)):
        ax1.plot(eval_pwm[i], colors[i], label=labels[i], linewidth=2.0)
        ax2.plot(eval_pwm[i], colors[i], label=labels[i], linewidth=2.0)
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF genome')
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
    ax2.get_yaxis().set_ticks([])

    # ax.set_xscale('log')
    # plt.xlim(0, len(eval_pwm[0]))

    # plt.savefig(dir_figures + parsed.figure_name + '_pwm.pdf', fmt='pdf')

    """ Bar plot at specific targets per tf level """
    fns = []

    # np original
    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_global_shrinkage/'
    dir_sub = 'analysis_compiled_chip_flynet_pwm/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_model_full.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_model_full.txt')
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_model_full_tf_merged_pid50.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_model_full_tf_merged_pid50.txt')

    dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_motif_incorporated/'

    # fire inference
    dir_sub = 'fire_motifs_np_tf_merged_dbd50_bin_20/analysis_compiled_chip_flynet_pwm/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_motif_net_fire_np_bin_20_resort_tf_merged.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_motif_net_fire_np_bin_20_resort_tf_merged.txt')

    dir_sub = 'fire_ortho_dmel+Dsim+Dsec_motifs_np_bin_20/analysis_compiled_chip_flynet_pwm/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')

    # cisbp inference
    dir_sub = 'cisbp_-2000_+200_fimo_dbd_cutoff_cellCycle_np_tf_merged_dbd50/analysis_compiled_chip_flynet_pwm/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_motif_net_dbd_cutoff_40_resort_tf_merged.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_motif_net_dbd_cutoff_40_resort_tf_merged.txt')

    # known motifs
    dir_sub = 'cisbp_-2000_+200_fimo_known_motif/analysis_compiled_chip_flynet_pwm/'
    fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_tf_merged_motif_net_known_motif_tf_merged_resort.txt')
    fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_tf_merged_motif_net_known_motif_tf_merged_resort.txt')

    # compute results
    [eval_chip, eval_pwm] = parse_binary_gold_standard(fns, parsed.eval_method)

    # plot figure
    xlevel = 100
    xlevel_index = xlevel/5-1 

    fig = plt.figure(num=None, figsize=(16,6), dpi=80)
    fig.subplots_adjust(wspace=.55)
    colors = ['0.5', 'k', 'b', 'r', 'm', 'c', 'g']
    labels = []
    labels.append('chance')
    labels.append('np')
    labels.append('np + tf_score_summing')
    labels.append('np + tf_score_summing\n + inferred_motif')
    labels.append('np + tf_score_summing\n + inferred_ortho_motif')
    labels.append('np + tf_score_summing\n + query_cisbp_motif')
    labels.append('np + tf_score_summing\n + known_motif')
    y_pos = numpy.arange(len(labels))

    plt.subplot(1,2,1)
    for i in y_pos:
        plt.barh(y_pos[i], eval_chip[:,xlevel_index][i], align='center', alpha=1, color=colors[i])
    plt.ylabel('Network mapping procedures')
    plt.xlabel('Interactions supported by ChIP (%)')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()

    ax = plt.subplot(1,2,2)
    for i in y_pos:
        plt.barh(y_pos[i], eval_pwm[:,xlevel_index][i], align='center', alpha=1, color=colors[i])
    plt.yticks(y_pos, labels)
    ax.set_yticklabels(labels, horizontalalignment='center', position=(-.25,.5))
    plt.xlabel('Interactions supported by PWM (%)')
    plt.gca().invert_yaxis()

    plt.savefig(dir_figures + parsed.figure_name + '_bar_' + str(xlevel) +'.pdf' , fmt='pdf')


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
