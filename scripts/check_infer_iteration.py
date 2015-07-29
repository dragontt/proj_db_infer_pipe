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

    # use only changed binned tf network scores for motif inference
    if parsed.iter_step == 'bin':
        counter = 0
        fns_curr = glob.glob(parsed.curr_iter + '*')
        for fn_curr in fns_curr:
            fn_prev = parsed.prev_iter + os.path.basename(fn_curr)
            if filecmp.cmp(fn_curr, fn_prev):
                os.system('rm '+ fn_curr)
                counter += 1
        print str(counter), 'unchanged TF removed.'

    # use only changed inferred motifs for FIMO scoring
    elif parsed.iter_step == 'motif':
        pass

    # link unchanged FIMO scores from prev. iteration for building motif net
    elif parsed.iter_step == 'fimo':
        pass 

    # incorrect argument
    else:
        sys.exit('Iteration step incorrect.')

def check_dir(directory):
    if not directory.endswith('/'):
        directory += '/'
    return directory

if __name__ == "__main__":
    main(sys.argv)
