#!/usr/bin/python

import sys
import argparse
import os.path
import glob
from scipy.stats import rankdata
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Compute rankings of fimo outputs")
    parser.add_argument('-i', '--input_dir', dest='input_dir', type=str)
    parser.add_argument('-o', '--output_dir', dest='output_dir', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    if not parsed.input_dir.endswith('/'):
        parsed.input_dir += '/'
    if not parsed.output_dir.endswith('/'):
        parsed.output_dir += '/'

    # get tf names
    tf_names = []
    filenames = glob.glob(parsed.input_dir + "*.summary")
    for fn in filenames:
        temp_name = os.path.basename(fn).split('.')
        tf_names.append(temp_name[0] + '.' + temp_name[1])

    # sort fimo rankings for each tf
    for i in range(len(tf_names)):
        lines = open(filenames[i], 'r').readlines()
        targets = []
	scores = []
        for line in lines:
            targets.append(line.split()[1].strip())
            scores.append(abs(float(line.split()[8].strip())))
        rankings = rankdata(scores)
	rankings = len(scores) + 1 - rankings
	sorted_indices = numpy.argsort(rankings)

        # write output file
        out_file = parsed.output_dir + tf_names[i]
        writer = open(out_file, 'w')
        for i in range(len(sorted_indices)):
            writer.write("%s\t%.2f\n" % (targets[sorted_indices[i]], rankings[sorted_indices[i]]))
        writer.close()

if __name__ == "__main__":
	main(sys.argv)
