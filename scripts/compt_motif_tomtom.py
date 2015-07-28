#!/usr/bin/python

"""
Align inferred motifs to known ones.
"""

import sys
import os
import argparse
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Align inferred motifs to known ones")
    parser.add_argument('-f', '-fn_known_motifs', dest='fn_known_motifs')
    parser.add_argument('-k', '-dir_known_pfms', dest='dir_known_pfms')
    parser.add_argument('-i', '-dir_inferred_pfms', dest='dir_inferred_pfms')
    parser.add_argument('-o', '-fn_output', dest='fn_output')
    parser.add_argument('-u', '-use_ortho', dest='use_ortho', type=bool, default=False)
    parser.add_argument('-s', '-species', dest='species')
    parser.add_argument('-t', '-fn_ortho', dest='fn_ortho')
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)
    parsed.dir_known_pfms = check_dir(parsed.dir_known_pfms)
    parsed.dir_inferred_pfms = check_dir(parsed.dir_inferred_pfms)

    # tomtom arguments
    bfile = '/home/mblab/ykang/proj_db_infer_pipe/resources/cisbp_all_species_bg_freq/Drosophila_melanogaster.bmf'
    E_value_thresh = 1

    # parse gene orthologs 
    if parsed.use_ortho:
        print 'True'
        ortho_dict = {}
        lines = open(parsed.fn_ortho, 'r').readlines()
        for i in range(5, len(lines)):
            line = lines[i].split('\t')
            if (line[6].startswith(parsed.species)):
                ortho_dict[line[0]] = line[5]

    # run tomtom 
    aligned_motifs = ['#Query ID\tTarget ID\tOptimal offset\tp-value\tE-value\tq-value\tOverlap\tQuery consensus\tTarget consensus\tOrientation\n']
    know_motifs = numpy.loadtxt(parsed.fn_known_motifs, dtype=str, skiprows=1)
    
    for i in range(len(know_motifs)):
        if parsed.use_ortho:
            if know_motifs[i,0] in ortho_dict.keys():
                tf = ortho_dict[know_motifs[i,0]]
        else:
            tf = know_motifs[i,0]

        # run tomtom on available inferred tfs
        if os.path.isfile(parsed.dir_inferred_pfms + tf):
            os.makedirs('tmp_tomtom')
            os.system('tomtom -no-ssc -oc tmp_tomtom -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh '+ str(E_value_thresh) +' -bfile '+ bfile +' '+ parsed.dir_inferred_pfms + tf +' '+ parsed.dir_known_pfms + know_motifs[i,1])
            if os.path.isfile('tmp_tomtom/'+ tf +'/tomtom.txt'):
                lines = open('tmp_tomtom/'+ tf +'/tomtom.txt', 'r').readlines()
                if len(lines) > 2:
                    aligned_motifs.append(lines[2])
            os.system('rm -r tmp_tomtom')

    # write tomtom summary
    writer = open(parsed.fn_output, 'w')
    for i in range(len(aligned_motifs)):
        writer.write('%s' % aligned_motifs[i])
    writer.close()

def check_dir(directory):
    if not directory.endswith('/'):
        directory += '/'
    return directory

if __name__ == "__main__":
    main(sys.argv)
