#!/usr/bin/python

"""
Prepare resources for Netprophet using external microarray epxression data. The output are: allowed.adj, conditions.environmental, data.expr, data.regulators.expr, data.pert, gids.fb, rids.fb, rderank.cuffdiff.adj
"""

import sys
import os
import argparse
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Prepare resources for network mapping")
    parser.add_argument('-p', '--ispert', dest='ispert', type=int)
    parser.add_argument('-e', '--fn_expr', dest='fn_expr', type=str)
    parser.add_argument('-g', '--fn_gids', dest='fn_gids', type=str)
    parser.add_argument('-r', '--fn_rids', dest='fn_rids', type=str)
    parser.add_argument('-o', '--dir_output', dest='dir_output', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    # load data
    gids = numpy.loadtxt(parsed.fn_gids, dtype=str)
    rids = numpy.loadtxt(parsed.fn_rids, dtype=str)
    expr_g = numpy.loadtxt(parsed.fn_expr, skiprows=1)
    conds = open(parsed.fn_expr, 'r').readline().strip().split('\t')

    # find the intersection of givens gids and rids
    rids = numpy.intersect1d(gids, rids)
    # size of network (N x M)
    N = len(rids)
    M = len(gids)

    # compute regulator expression data, and compute allowed matrix (disallowing self-regulation)
    expr_r = numpy.zeros((N,len(conds)))
    allowed = numpy.ones((N,M))
    for i in range(N):
        j = numpy.where(gids == rids[i])[0][0]
        expr_r[i,:] = expr_g[j,:]
        allowed[i,j] = 0

    # create empty perturbation and differential expreesion matrices if not gene perturbation occured
    pert = numpy.zeros(expr_g.shape)
    if parsed.ispert == 0:
        de = numpy.zeros((N,M))
    else:
        for j in range(len(conds)):
            i = numpy.where(gids == conds[j])[0]
            if len(i) != 0:
                pert[i[0],j] = 1

    # check output directory
    if not os.path.exists(parsed.dir_output):
        os.makedirs(parsed.dir_output)
    if not parsed.dir_output.endswith("/"):
        parsed.dir_output += "/"

    # write data: allowed.adj, conditions.environmental, data.expr, data.regulators.expr, data.pert, gids.fb, rids.fb, rderank.cuffdiff.adj
    numpy.savetxt(parsed.dir_output+"allowed.adjmtr", allowed, fmt='%d')
    numpy.savetxt(parsed.dir_output+"conditions.environmental", conds, fmt='%s')
    write_matrix(parsed.dir_output+"data.expr", expr_g)
    write_matrix(parsed.dir_output+"data.regulators.expr", expr_r)
    numpy.savetxt(parsed.dir_output+"data.pert", pert, fmt='%d')
    numpy.savetxt(parsed.dir_output+"gids.fb", gids, fmt='%s')
    numpy.savetxt(parsed.dir_output+"rids.fb", rids, fmt='%s')
    if parsed.ispert == 0:
        numpy.savetxt(parsed.dir_output+"rderank.cuffdiff.adj", de, fmt='%d')
    else:
        print 'DE adjmtr is not computed; prepare DE information separately.'

def write_matrix(fn, mtr):
    writer = open(fn, 'w')
    for i in range(mtr.shape[0]):
        for j in range(mtr.shape[1]):
            if mtr[i,j] == 0:
                writer.write("0\t")
            else:
                writer.write("%.10f\t" % mtr[i,j])
        writer.write("\n")
    writer.close()

if __name__ == "__main__":
    main(sys.argv)
