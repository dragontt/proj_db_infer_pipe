#!/usr/bin/python
""" 
Combine mutliple networks using arithmetic mean, with weights assgined to the networks.
Input: the first (N-1) files of networks (adjmtr)
Output: the last one
"""

import sys
import argparse
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine network using geometric mean")
    parser.add_argument("ARGS", help="Args to store as Gist", nargs="+")
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)
    
    # get input networks
    args = parsed.ARGS
    networks = []
    for i in range((len(args)-1)/2):
        x = numpy.abs(numpy.loadtxt(args[i]))
        networks.append(x)
    weights = []
    for i in range((len(args)-1)/2, len(args)-1):
        weights.append(float(args[i]))
    combined = numpy.average(networks, weights=weights, axis=0)

    # write combined network
    fn_output = parsed.ARGS[len(args)-1]
    # numpy.savetxt(fn_output, combined, fmt="%.10f", delimiter="\t", newline="\n")
    write_adjmtr(fn_output, combined)

def write_adjmtr(fn, adjmtr):
    writer = open(fn, "w")
    for i in range(len(adjmtr)):
        for j in range(len(adjmtr[i])):
            if adjmtr[i,j] == 0:
                writer.write("0\t")
            else:
                writer.write("%0.10f\t" % adjmtr[i,j])
        writer.write("\n")
    writer.close()

if __name__ == "__main__":
    main(sys.argv)