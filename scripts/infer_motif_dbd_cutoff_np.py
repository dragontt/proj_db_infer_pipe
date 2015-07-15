#!/usr/bin/python

"""
Infer a motif by first filtering database motifs based on the DNA binding domain percent identity (pid) within certain pid cutoff, and obtain the motif with highest rank order correlation of FIMO scan scores and the inferred network regulator-target scores.
"""

import sys
import argparse
import numpy
import glob 
import os.path

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Infer a motif using both DBD pid of k-NN and highest correlation with network edge score.")
    parser.add_argument('-a', '--dir_dbd_pid', dest='dir_dbd_pid', type=str)
    parser.add_argument('-c', '--dir_rk_corr', dest='dir_rk_corr', type=str)
    parser.add_argument('-p', '--param_cutoff', dest='param_cutoff', type=float, default=40)
    parser.add_argument('-o', '--dir_output', dest='dir_output', type=str)
    parser.add_argument('-d', '--fn_conv_dbd2rid', dest='fn_conv_dbd2rid', type=str, default='/home/mblab/ykang/proj_db_infer_pipe/resources/fly_aa_seq/pids.dbd2fbgn')
    parser.add_argument('-l', '--fn_not_use', dest='fn_not_use', type=str, default='/home/mblab/ykang/proj_db_infer_pipe/resources/cisbp_1.01/cisbp_motifs_Drosophila_melanogaster.txt')
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)
    
    # prepare data
    parsed.dir_dbd_pid = check_dir(parsed.dir_dbd_pid)
    parsed.dir_rk_corr = check_dir(parsed.dir_rk_corr)
    parsed.dir_output = check_dir(parsed.dir_output)

    dbds = get_basenames(parsed.dir_dbd_pid)
    motifs_not_use = parse_list(parsed.fn_not_use)
    dbd2fbgn = parse_dict(parsed.fn_conv_dbd2rid, 'str')

    regulators = list(set(dbd2fbgn.values())) 
    inferred_motif = {}
    inferred_pid = {}
    inferred_corr = {}

    # infer the motif using dbd and np corr
    for dbd in dbds:
        regulator = dbd2fbgn[dbd]
        fn_regulator = parsed.dir_rk_corr + regulator

        if os.path.isfile(fn_regulator):
            # sys.stdout.write("%s, " % regulator)
            reg_corrs = parse_dict(fn_regulator, 'float')

            candidate_motif = []
            candidate_pid = []
            candidate_corr = []

            lines = open(parsed.dir_dbd_pid + dbd, 'r').readlines()
            for i in range(len(lines)):
                line = lines[i].strip().split()
                motif = line[0].split(':')[0]
                pid = float(line[1])

                if (pid >= parsed.param_cutoff):
                    if (not motif in motifs_not_use) and (motif in reg_corrs.keys()):
                        candidate_motif.append(motif)
                        candidate_pid.append(pid)
                        candidate_corr.append(reg_corrs[motif])
                else: 
                    break

            if len(candidate_corr) > 0:
                index_max = numpy.argmax(candidate_corr)

                if not regulator in inferred_motif.keys():
                    inferred_motif[regulator] = candidate_motif[index_max]
                    inferred_pid[regulator] = candidate_pid[index_max]
                    inferred_corr[regulator] = candidate_corr[index_max]
                else:
                    if candidate_corr[index_max] > inferred_corr[regulator]:
                        inferred_motif[regulator] = candidate_motif[index_max]
                        inferred_pid[regulator] = candidate_pid[index_max]
                        inferred_corr[regulator] = candidate_corr[index_max]

    # write data
    writer = open(parsed.dir_output + "motifs_dbd_cutoff_" + str(parsed.param_cutoff) + ".txt", "w")
    writer.write("#query\tinferred_motif(s)\tdbd_pid\tnp_corr\n")
    for regulator in sorted(inferred_motif.keys()):
        writer.write("%s\t" % regulator)
        if len(inferred_motif[regulator]) == 0:
            writer.write("None\tNone\tNone\n")
        else:
            writer.write("%s\t%.5f\t%.5f\n" % (inferred_motif[regulator], inferred_pid[regulator], inferred_corr[regulator]))
    writer.close()

def check_dir(file_dir):
    if not file_dir.endswith('/'):
        file_dir += '/'
    return file_dir

def get_basenames(file_dir):
    names = []
    filenames = glob.glob(file_dir + "*")
    for filename in filenames:
	filename = os.path.basename(filename)
        if not filename.startswith('_'):
            names.append(filename)
    return names

def parse_dict(fn, dict_type):
    out = {}
    lines = open(fn, 'r').readlines()
    if dict_type == 'str':
        for line in lines:
            temp = line.strip().split()
            out[temp[0]] = temp[1]
    elif dict_type == 'float':
        for line in lines:
            temp = line.strip().split()
            if temp[1] != 'nan':
                out[temp[0]] = float(temp[1])
    return out

def parse_list(fn):
    names = set()
    lines = open(fn, 'r').readlines()
    for line in lines:
        names.add(line.strip())
    return names

if __name__ == "__main__":
    main(sys.argv)
