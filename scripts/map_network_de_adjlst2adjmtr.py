#!/usr/bin/python
import sys
import argparse
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Map DE adjlst to adjmtr using -log(pval).")
    parser.add_argument('-l', '--fn_adjlst', dest='fn_adjlst', type=str)
    parser.add_argument('-m', '--fn_adjmtr', dest='fn_adjmtr', type=str)
    parser.add_argument('-r', '--fn_rids', dest='fn_rids', type=str)
    parser.add_argument('-g', '--fn_gids', dest='fn_gids', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)
    convert_adjlst2adjmtr(parsed.fn_adjlst, parsed.fn_adjmtr, parsed.fn_rids, parsed.fn_gids, 1)

def convert_adjlst2adjmtr(fn_adjlst, fn_adjmtr, fn_rids, fn_gids, skiprows):
    rids = numpy.loadtxt(fn_rids, dtype=str)
    gids = numpy.loadtxt(fn_gids, dtype=str)
    lines = open(fn_adjlst, 'r').readlines()
    adjmtr = numpy.zeros((len(rids), len(gids)))   
    for i in range(skiprows,len(lines)):
        line = lines[i].strip().split()
        de_score = - numpy.log2(float(line[4]))
        r_index = numpy.where(rids == line[0])[0]
        g_index = numpy.where(gids == line[1])[0]
        if len(r_index) != 0 and len(g_index) != 0:
            adjmtr[r_index[0], g_index[0]] = de_score

    write_matrix(fn_adjmtr, adjmtr)   

def parse_id_index(lst, ids):
    out_dict = {}
    for name in lst:
        index = numpy.where(ids == name)[0]
        if len(index) != 0:
            out_dict[name] = index[0]
    return out_dict

def write_matrix(fn, mtr):
    writer = open(fn, 'w')
    for i in range(mtr.shape[0]):
        for j in range(mtr.shape[1]):
            if mtr[i,j] == 0:
                writer.write("0\t")
            else:
                writer.write("%.15f\t" % mtr[i,j])
        writer.write("\n")
    writer.close()

if __name__ == "__main__":
    main(sys.argv)