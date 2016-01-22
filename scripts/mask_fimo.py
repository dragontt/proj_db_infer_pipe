#!/usr/bin/python
"""
Mask FIMO scan output by the conserved element loci masked promoters.
"""

import sys
import os
import argparse
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Mask FIMO scan output by the conserved element masked promoters.")
    parser.add_argument('-i', '-fn_input_fimo', dest='fn_input_fimo')
    parser.add_argument('-a', '-fn_mask', dest='fn_mask')
    parser.add_argument('-m', '-mask_type', dest='mask_type', default='mask1')
    parser.add_argument('-o', '-fn_output_fimo', dest='fn_output_fimo')
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)

    # read files
    prom_mask = {}
    lines = open(parsed.fn_mask, 'r').readlines()
    for i in range(len(lines)):
        line = lines[i].split('\t')
        gid = line[0]
        prom_mask[gid] = []
        if len(line) > 1:
            for j in range(len(line)-1):
                prom_mask[gid].append(numpy.array(line[j+1].split(' '), dtype=int))
        prom_mask[gid] = numpy.array(prom_mask[gid])

    # filter fimo file
    lines = open(parsed.fn_input_fimo, 'r').readlines()

    writer = open(parsed.fn_output_fimo, 'w')
    writer.write('%s' % lines[0])
    
    for i in range(1, len(lines)):
        line = lines[i].split()
        (gid, start, stop) = (line[1], int(line[2]), int(line[3]))

        if (gid in prom_mask.keys()) and (len(prom_mask[gid]) > 0):
            
            if parsed.mask_type ==  'mask1':
                # hit within the masked loci
                inds = numpy.intersect1d(numpy.where(prom_mask[gid][:,0]<=start)[0],
                    numpy.where(prom_mask[gid][:,1]>=stop)[0])

            elif parsed.mask_type == 'mask2':
                # hit intersects with masked loci
                inds = numpy.union1d(
                    numpy.intersect1d(numpy.where(prom_mask[gid][:,0]<=start)[0],
                    numpy.where(prom_mask[gid][:,1]>=start)[0]),
                    numpy.intersect1d(numpy.where(prom_mask[gid][:,0]<=stop)[0],
                    numpy.where(prom_mask[gid][:,1]>=stop)[0]))

            elif parsed.mask_type == 'mask3':
                # hit within certain shift of the masked loci 
                shift = 27 # 3x of the motif width
                inds = numpy.intersect1d(numpy.where(prom_mask[gid][:,0]-shift<=start)[0],
                    numpy.where(prom_mask[gid][:,1]+shift>=stop)[0])

            else:
                sys.stderr.write('Unspecified mask type')

            if len(inds) > 0:
                writer.write('%s' % lines[i])

    writer.close()

if __name__ == "__main__":
    main(sys.argv)
