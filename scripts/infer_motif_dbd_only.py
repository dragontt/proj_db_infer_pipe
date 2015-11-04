#!/usr/bin/python

"""
Infer a motif based on the highest DNA binding domain percent identity (pid), with pid bottom line of 40%.
"""

import sys
import argparse
import numpy
import glob 
import os.path

PID_THRESHOLD = 40

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Infer a motif only using dbd pid.")
    parser.add_argument('-a', '--dir_dbd_pid', dest='dir_dbd_pid', type=str)
    parser.add_argument('-o', '--dir_output', dest='dir_output', type=str)
    parser.add_argument('-u', '--use_conv', dest='use_conv', type=bool, default=False)
    parser.add_argument('-d', '--fn_conv', dest='fn_conv', type=str, default='/home/mblab/ykang/proj_db_infer_pipe/resources/fly_aa_seq/pids.dbd2fbgn')
    parser.add_argument('-l', '--fn_not_use', dest='fn_not_use', type=str, default='/home/mblab/ykang/proj_db_infer_pipe/resources/cisbp_1.01/cisbp_motifs_Drosophila_melanogaster.txt')
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)
    
    # prepare data
    parsed.dir_dbd_pid = check_dir(parsed.dir_dbd_pid)
    parsed.dir_output = check_dir(parsed.dir_output)

    dbds = get_basenames(parsed.dir_dbd_pid)
    motifs_not_use = parse_list(parsed.fn_not_use)	
    if parsed.use_conv:
        dbd2fbgn = parse_dict(parsed.fn_conv, 'str')
        regulators = list(set(dbd2fbgn.values())) 
    else:
        regulators = dbds
        dbd2fbgn = {}
        for regulator in regulators:
            dbd2fbgn[regulator] = regulator
 
    inferred_motif = {}
    inferred_pid = {}

    # infer the motif with highest dbd pid and inlcude the tied ones
    for dbd in dbds:
        regulator = dbd2fbgn[dbd]
        lines = open(parsed.dir_dbd_pid + dbd, 'r').readlines()        
        for i in range(len(lines)):
            line = lines[i].strip().split()
            motif = line[0].split(':')[0]
            pid = float(line[1])
            if (not motif in motifs_not_use) and pid > PID_THRESHOLD:
                if not regulator in inferred_motif.keys():
                    inferred_motif[regulator] = [motif]
                    inferred_pid[regulator] = pid
                else:
                    if pid == inferred_pid[regulator]:
                        inferred_motif[regulator].append(motif)
                    elif pid > inferred_pid[regulator]:
                        inferred_motif[regulator] = [motif]
                        inferred_pid[regulator] = pid
                    else:
                        break
    for key in inferred_motif.keys():
	   inferred_motif[key] = list(set(inferred_motif[key]))

    # write data
    writer = open(parsed.dir_output + "motifs_dbd_only.txt" , "w")
    writer.write("#query\tinferred_motif(s)\tdbd_pid\n")
    for regulator in sorted(inferred_motif.keys()):
        writer.write("%s\t" % regulator)
        inferred_length = len(inferred_motif[regulator])
        if inferred_length == 0:
            writer.write("None\tNone\n")
        elif inferred_length == 1:
            writer.write("%s\t%.5f\n" % (inferred_motif[regulator][0], inferred_pid[regulator]))
        else:
            for i in range(inferred_length-1):
                writer.write("%s, " % inferred_motif[regulator][i])
            writer.write("%s\t%.5f\n" % (inferred_motif[regulator][inferred_length-1], inferred_pid[regulator]))
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
