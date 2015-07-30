#!/usr/bin/python

"""
Check motif and network inference iteration. The steps to check as duplicate comparison to the previous iteration includes:
    1. bin:     use only changed binned tf network scores for motif inference
    2. motif:   use only changed inferred motifs for FIMO scoring
    3. fimo:    link unchanged FIMO scores from prev. iteration for building motif net 
"""

import sys
import os
import glob
import argparse
import numpy
import filecmp

iter_steps = ['bin', 'motif', 'fimo']

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Align inferred motifs to known ones")
    parser.add_argument('-s', '-iter_step', dest='iter_step', help='options: %s' % iter_steps)
    parser.add_argument('-pi', '-prev_iter', dest='prev_iter')
    parser.add_argument('-ci', '-curr_iter', dest='curr_iter')
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)
    if parsed.iter_step == 'bin' or parsed.iter_step == 'fimo':
        parsed.prev_iter = check_dir(parsed.prev_iter)
        parsed.curr_iter = check_dir(parsed.curr_iter)

    counter = 0

    # use only changed binned tf network scores for motif inference
    if parsed.iter_step == 'bin':
        fns_curr = glob.glob(parsed.curr_iter + '*')
        for fn_curr in fns_curr:
            fn_prev = parsed.prev_iter + os.path.basename(fn_curr)
            if filecmp.cmp(fn_curr, fn_prev):
                os.system('rm '+ fn_curr)
                counter += 1
        print str(counter), 'unchanged TF removed.'

    # use only changed inferred motifs for FIMO scoring
    elif parsed.iter_step == 'motif':
        motifs_dict = {}
        motifs_prev = numpy.loadtxt(parsed.prev_iter, dtype=str)
        for i in range(len(motifs_prev)):
            motifs_dict[motifs_prev[i,0]] = motifs_prev[i,1]
        
        motifs_out = numpy.empty([0,2], dtype=str)
        motifs_curr = numpy.loadtxt(parsed.curr_iter, dtype=str)
        for i in range(len(motifs_curr)):
            if (motifs_curr[i,0] in motifs_dict.keys()) and (motifs_curr[i,1] != motifs_dict[motifs_curr[i,0]]):
                motifs_out = numpy.vstack((motifs_out, motifs_curr[i,:]))
        numpy.savetxt(parsed.curr_iter, motifs_out, fmt='%s')
        print str(counter), 'unchanged TF removed.'

    # link unchanged FIMO scores from prev. iteration for building motif net
    elif parsed.iter_step == 'fimo':
        motifs_set = set()
        fns_curr = glob.glob(parsed.curr_iter + '*')
        for fn in fns_curr:
            motifs_set.add(os.path.basename(fn))
        
        fns_prev = glob.glob(parsed.prev_iter + '*')
        for fn in fns_prev:
            fn_basename = os.path.basename(fn)
            if not fn_basename in motifs_set:
                os.system('ln -s '+ os.path.abspath(parsed.prev_iter) +'/'+ fn_basename +' '+ parsed.curr_iter + fn_basename)
                counter += 1
        print str(counter), 'prev avail TF linked.'

    # incorrect argument
    else:
        sys.exit('Iteration step incorrect.')

def check_dir(directory):
    if not directory.endswith('/'):
        directory += '/'
    return directory

if __name__ == "__main__":
    main(sys.argv)
