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
    parser.add_argument('-n', '-figure_name', dest='figure_name', default='temp_result')
    parser.add_argument('-r', '-range', dest='range', default='top4to40k')
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

    dir_figures = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_analysis_results/'
    figure_title = 'Fly CellCycle Network: Combined with Inferred Motifs'

    if parsed.eval_type == "flynet":

        # evaluate chip and pwm supports on binary gold standard 
        # file initialization
        fns = []

        # np original
        dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_global_shrinkage/'
        dir_sub = 'analysis_compiled_chip_flynet_pwm/'
        fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_model_full.txt')
        fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_model_full.txt')
        fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_model_full_tf_merged_pid50.txt')
        fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_model_full_tf_merged_pid50.txt')

        # known motifs
        dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/fly_network_cellCycle_motif_incorporated/'
        dir_sub = 'cisbp_-2000_+200_fimo_known_motif/analysis_compiled_chip_flynet_pwm/'
        # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_motif_net_known_motif_resort.txt')
        # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_motif_net_known_motif_resort.txt')
        # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_tf_merged_motif_net_known_motif_tf_merged_resort.txt')
        # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_tf_merged_motif_net_known_motif_tf_merged_resort.txt')

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
        fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')
        fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')
        dir_sub = 'fire_ortho_dmel+Dsim+Dsec+Dyak+Dere+Dana_motifs_np_bin_20/analysis_compiled_chip_flynet_pwm/'
        fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')
        fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_ortho_motif_net_tf_merged_resort.txt')

        # cisbp inference
        dir_sub = 'cisbp_-2000_+200_fimo_dbd_cutoff_cellCycle_np/analysis_compiled_chip_flynet_pwm/'
        # fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_motif_net_dbd_cutoff_40.0_resort_0.txt')
        # fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_motif_net_dbd_cutoff_40.0_resort_0.txt')

        dir_sub = 'cisbp_-2000_+200_fimo_dbd_cutoff_cellCycle_np_tf_merged_dbd50/analysis_compiled_chip_flynet_pwm/'
        fns.append(dir_network + dir_sub + 'analysis_chip_support.'+ parsed.range +'.combined_network_np_motif_net_dbd_cutoff_40_resort_tf_merged.txt')
        fns.append(dir_network + dir_sub + 'analysis_pwm_support.'+ parsed.range +'.combined_network_np_motif_net_dbd_cutoff_40_resort_tf_merged.txt')

        # figure setup
        # colors = ['k:', 'k', 'k--', 'm', 'm--', 'r', 'r--', 'g', 'g--', 'g:', 'b', 'b--' ]
        colors = ['k:', 'k--', 'k', 'r', 'g', 'c', 'b']
        labels = []
        labels.append('chance')
        labels.append('np')
        labels.append('np_tf_merged')
        # labels.append('np + known_motifs')
        # labels.append('np_tf_merged + known_motifs')
        # labels.append('np + fire_dmel_motifs')
        labels.append('np_tf_merged + fire_dmel_motifs')
        # labels.append('np_tf_merged + fire_indiv_ortho_motifs')
        labels.append('np_tf_merged + fire_comb_ortho_motifs _2species')
        labels.append('np_tf_merged + fire_comb_ortho_motifs _5species')
        # labels.append('np + cisbp_motifs')
        labels.append('np_tf_merged + cisbp_motifs')

        if parsed.range == 'top1to10k':
            x_ticks = ['1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k']
        elif parsed.range == 'top4to40k':
            x_ticks = ['4k', '8k', '12k', '16k', '20k', '24k', '28k', '32k', '36k', '40k']
        elif parsed.range == 'top12to120k':
            x_ticks = ['12k', '24k', '36k', '48k', '60k', '72k', '84k', '96k', '108k', '120k']
        elif parsed.range == 'top80to800k':
            x_ticks = ['80k', '160k', '240k', '320k', '400k', '480k', '560k', '640k', '720k', '800k']
        elif parsed.range == 'top120to1200k':
            x_ticks = ['120k', '240k', '360k', '480k', '600k', '720k', '840k', '960k', '1080k', '1200k']
        else:
            x_ticks = []

        # compute chip and pwm supports
        [eval_chip, eval_pwm] = parse_binary_gold_standard(fns, parsed.eval_method)

    elif parsed.eval_type == "binding_indep":
        pass

    # plot figures
    # fig = plt.figure(num=None, figsize=(25,8), dpi=80)
    # fig.subplots_adjust(right=.7)

    fig = plt.figure(num=None, figsize=(18,8), dpi=80)

    plt.subplot(1,2,1)
    for i in range(len(eval_chip)):
        plt.plot(eval_chip[i], colors[i], label=labels[i])
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by ChIP')
    plt.xlim(-1, len(eval_chip[0]))
    if parsed.range == 'top1to10k':
        plt.ylim(0, 1)
    elif parsed.range == 'top4to40k':
        plt.ylim(0, .6)
    else:
        plt.ylim(0, .6)

    plt.subplot(1,2,2)
    for i in range(len(eval_pwm)):
        plt.plot(eval_pwm[i], colors[i], label=labels[i])
    plt.xticks(range(len(eval_pwm[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by PWM')
    plt.xlim(-1, len(eval_pwm[0]))
    plt.ylim(0, .2)
    plt.legend(loc="upper right")
    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    plt.suptitle(figure_title)

    plt.savefig(dir_figures + parsed.figure_name + '.pdf', fmt='pdf')
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
