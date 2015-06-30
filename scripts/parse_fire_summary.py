#!/usr/bin/python

import sys
import argparse
import glob
import os.path
import os

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Summarize inferred FIRE motifs.")
    parser.add_argument('-i', '--dir_input', dest='dir_input', type=str)
    parser.add_argument('-o', '--fn_output', dest='fn_output', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    # check directory
    if not parsed.dir_input.endswith("/"):
        parsed.dir_input += "/"

    writer = open(parsed.fn_output, "w")

    # parse the lists of tfs and inferred motifs
    fns = glob.glob(parsed.dir_input + "*_FIRE")
    for fn in fns:
        tf = os.path.basename(fn).split('_')[0]
        # parse FIRE summary file  
        fsum = fn + "/DNA/" + tf + ".summary"
        if os.path.isfile(fsum) and os.stat(fsum).st_size > 0:
            lines = open(fsum, "r").readlines()
            motif = lines[0].split()[0]
            writer.write("%s\t%s\n" % (tf, motif))
    writer.close()

if __name__ == "__main__":
    main(sys.argv)
