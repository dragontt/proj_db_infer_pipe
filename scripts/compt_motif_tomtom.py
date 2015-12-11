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
    parser.add_argument('-k', '-fn_known_pfms', dest='fn_known_pfms')
    parser.add_argument('-i', '-dir_inferred_pfms', dest='dir_inferred_pfms')
    parser.add_argument('-e', '-E_value_thresh', dest='E_value_thresh')
    parser.add_argument('-b', '-fn_background_freq', dest='fn_background_freq')
    parser.add_argument('-o', '-fn_output', dest='fn_output')
    parser.add_argument('-u', '-use_ortho', dest='use_ortho', type=bool, default=False)
    parser.add_argument('-s', '-species', dest='species', default='Dmel')
    parser.add_argument('-t', '-fn_ortho', dest='fn_ortho')
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)
    parsed.dir_inferred_pfms = check_dir(parsed.dir_inferred_pfms)

    # tomtom arguments
    # bfile = '/home/mblab/ykang/proj_db_infer_pipe/resources/cisbp_all_species_bg_freq/Drosophila_melanogaster.bmf'
    if parsed.fn_background_freq == None:
        sys.exit("No background frequency file loaded.")
    else:
        bfile = parsed.fn_background_freq
    # E_value_thresh = 1

    # parse gene orthologs
    if parsed.use_ortho:
        ortho_dict = {}
        lines = open(parsed.fn_ortho, 'r').readlines()
        for i in range(5, len(lines)):
            line = lines[i].split('\t')
            if (line[6].startswith(parsed.species)):
                ortho_dict[line[0]] = line[5]

    # run tomtom 
    results = numpy.array([['#TF','ortho_TF','Query_ID','Target_ID','Optimal_offset','p-value','E-value','q-value','Overlap','Query_consensus','Target_consensus','Orientation']])
    known_motifs = numpy.loadtxt(parsed.fn_known_motifs, dtype=str, skiprows=1)
    
    os.system('mkdir -p tomtom_output_Eval_' + parsed.E_value_thresh)

    for i in range(len(known_motifs)):
        if parsed.use_ortho:
            if known_motifs[i,0] in ortho_dict.keys():
                tf = ortho_dict[known_motifs[i,0]]
        else:
            tf = known_motifs[i,0]

        # run tomtom on available inferred tfs
        if os.path.isfile(parsed.dir_inferred_pfms + tf):
            tmp_dir_tomtom = 'tomtom_output/'+ parsed.species +'_'+ tf
            os.makedirs(tmp_dir_tomtom)
            os.system('tomtom -no-ssc -oc '+ tmp_dir_tomtom +' -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh '+ str(parsed.E_value_thresh) +' -bfile '+ bfile +' '+ parsed.dir_inferred_pfms + tf +' '+ parsed.fn_known_pfms)
            if os.path.isfile(tmp_dir_tomtom +'/tomtom.txt'):
                tomtom_result = numpy.loadtxt(tmp_dir_tomtom +'/tomtom.txt', dtype=str, skiprows=1, delimiter='\t', ndmin=2)
                # search for the known motif for that tf
                if (tomtom_result.size > 0) and (known_motifs[i,1] in tomtom_result[:,1]):
                    index = numpy.where(tomtom_result[:,1] == known_motifs[i,1])[0]
                    if index.size > 0:
                        results = numpy.vstack((results, numpy.append([known_motifs[i,0], tf], tomtom_result[index[0],:])))
                # os.system('echo Dmel: '+ known_motifs[i,0] +' cisbp: '+ known_motifs[i,1] +' ortho: '+ tf + ' >> log')
                # os.system('cat '+ tmp_dir_tomtom +'/tomtom.txt >> log')

    # os.system('rm -r tomtom_output')

    # write results
    numpy.savetxt(parsed.fn_output, results, fmt='%s', delimiter='\t')

def check_dir(directory):
    if not directory.endswith('/'):
        directory += '/'
    return directory

if __name__ == "__main__":
    main(sys.argv)
