#!/usr/bin/python

"""
Query a list of database inference results, combine the inferred pwm fimo scan scores,
which were preprocessed in the motif inference step.
Inference methods: (1) cisbp, (2) fire
"""

import sys
import argparse
import os.path
import numpy
from scipy.stats.mstats import gmean

inference_methods = ['cisbp', 'fire']

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine inferred pwm fimo scan scores.")
    parser.add_argument('-m', '--inference_method', dest='inference_method', type=str, default='cisbp', help='options: %s' % inference_methods)
    parser.add_argument('-i', '--fn_infer', dest='fn_infer', type=str)
    parser.add_argument('-r', '--fn_rids', dest='fn_rids', type=str)
    parser.add_argument('-g', '--fn_gids', dest='fn_gids', type=str)
    parser.add_argument('-f', '--dir_fimo', dest='dir_fimo', type=str)    
    parser.add_argument('-o', '--fn_adjmtr', dest='fn_adjmtr', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    # check directory
    if not parsed.dir_fimo.endswith("/"):
        parsed.dir_fimo += "/"

    # get the lists of tfs, targets
    rids = numpy.loadtxt(parsed.fn_rids, dtype=str)
    gids = numpy.loadtxt(parsed.fn_gids, dtype=str)
    adjmtr = numpy.zeros([len(rids), len(gids)])

    # build the adjmtr
    lines = open(parsed.fn_infer, "r").readlines()
    for i in range(1, len(lines)):

        # get inferred tf and database motif
        linesplit = lines[i].strip().split('\t')
        infer_tf = linesplit[0]
        if parsed.inference_method == 'cisbp':
            infer_motifs = linesplit[1].split(',')
        elif parsed.inference_method == 'fire':
            infer_motifs = linesplit[0].split(',')
        else:
            sys.exit("Inference method not specified.")

        # get fimo scores for the inferred motif
        index = numpy.where(rids == infer_tf)[0]
        if len(index) > 0:
            if len(infer_motifs) > 1:
                temp_mtr = numpy.zeros([len(infer_motifs), len(gids)])
                for j in range(len(infer_motifs)):
                    fn_motif = parsed.dir_fimo + infer_motifs[j] + ".summary"
                    if os.path.isfile(fn_motif):
                        dict_scores = get_fimo_scores(fn_motif)
                        for k in range(len(gids)):
                            t = gids[k]
                            temp_mtr[j, k] = dict_scores[t] if t in dict_scores else 0
                adjmtr[index[0], :] = gmean(temp_mtr).data

            else:
                fn_motif = parsed.dir_fimo + infer_motifs[0] + ".summary"   
                if os.path.isfile(fn_motif):
                    dict_scores = get_fimo_scores(fn_motif)
                    for j in range(len(gids)):
                        t = gids[j]
                        adjmtr[index[0], j] = dict_scores[t] if t in dict_scores else 0
        
    # write adjmtr file
    write_adjmtr(adjmtr, parsed.fn_adjmtr)

def get_fimo_scores(fn):
    lines = open(fn, "r").readlines()
    d = {}
    for line in lines:
        temp_linesplit = line.strip().split()
        temp_name = temp_linesplit[1]
        # temp_score = (float(temp_linesplit[3]) + float(temp_linesplit[5]))/2
        temp_score = max([float(temp_linesplit[3]), float(temp_linesplit[5])])
        d[temp_name] = temp_score
    return d

def write_adjmtr(adjmtr, fn):
    writer = open(fn, "w")
    for i in range(len(adjmtr)):
        for j in range(len(adjmtr[i])):
            if adjmtr[i,j] == 0:
                writer.write("0.\t")
            else:
                writer.write("%0.10f\t" % adjmtr[i,j])
        writer.write("\n")
    writer.close()

if __name__ == "__main__":
    main(sys.argv)
