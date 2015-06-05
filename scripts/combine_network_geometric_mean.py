#!/usr/bin/python
""" 
Combine mutliple networks using geometric mean. All edge scores are converted to absolute value and zero score is replaced with the minimal absolute edge score of each network.
Input: the first (N-1) files of networks (adjmtr)
Output: the last one
"""

import sys
import argparse
import numpy
from scipy.stats.mstats import gmean

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine network using geometric mean")
    parser.add_argument("FILE", help="File to store as Gist", nargs="+")
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)
    
    # get input networks
    fns = parsed.FILE
    networks = []
    for i in range(len(fns)-1):
        x = numpy.abs(numpy.loadtxt(fns[i]))
        x += numpy.min(x[numpy.nonzero(x)])
        networks.append(x)
    combined = gmean(networks)

    # write combined network
    fn_output = parsed.FILE[len(parsed.FILE)-1]
    numpy.savetxt(fn_output, combined, fmt="%.10f", delimiter="\t", newline="\n")

if __name__ == "__main__":
    main(sys.argv)