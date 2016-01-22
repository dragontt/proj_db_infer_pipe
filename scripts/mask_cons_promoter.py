#!/usr/bin/python
"""
Mask promoters using UCSC conserved elements.
"""

import sys
import os
import argparse
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Mask promoters using UCSC conserved elements.")
    parser.add_argument('-g', '-fn_gids', dest='fn_gids')
    parser.add_argument('-a', '-fn_gene_annot', dest='fn_gene_annot')
    parser.add_argument('-e', '-fn_conselem', dest='fn_conselem')
    parser.add_argument('-l', '-promoter_length', dest='promoter_length', type=int)
    parser.add_argument('-o', '-fn_output_mask', dest='fn_output_mask')
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)

    # read files
    gids = numpy.loadtxt(parsed.fn_gids, dtype=str)
    gene_annot = numpy.loadtxt(parsed.fn_gene_annot, dtype=str, usecols=(1,2,6), delimiter='\t')
    cons_elem = numpy.loadtxt(parsed.fn_conselem, dtype=str, usecols=(1,2,3), delimiter='\t')
    prom_len = parsed.promoter_length
    
    # mask promoters 
    prom_mask = {}
    chrms = numpy.unique(gene_annot[:,1])
    for chrm in chrms:
        subset_gene_annot = gene_annot[numpy.where(gene_annot[:,1]==chrm)[0]][:,(0,2)]
        subset_cons_elem = cons_elem[numpy.where(cons_elem[:,0]==chrm)[0]][:,(1,2)]
        subset_cons_elem = numpy.array(subset_cons_elem, dtype=int)
        
        # mask promoter of each target gene
        for i in range(len(subset_gene_annot)):
            (gid, tss) = (subset_gene_annot[i,0], int(subset_gene_annot[i,1]))
            
            if gid in gids:
                (prom_start, prom_end) = (tss-prom_len, tss-1)
                inds = numpy.intersect1d(numpy.where(subset_cons_elem[:,0]>=prom_start)[0],
                    numpy.where(subset_cons_elem[:,1]<=prom_end)[0])
                # masked loci are the relative positions within promoter
                prom_mask[gid] = subset_cons_elem[inds] - tss + prom_len

    # write the mask
    writer = open(parsed.fn_output_mask, 'w')
    for gid in gids: 
        writer.write('%s' % gid)
        if gid in prom_mask.keys():
            loci = prom_mask[gid]
            if len(loci)>0:
                for i in range(len(loci)):
                    writer.write('\t%d %d' % (loci[i,0], loci[i,1]))
            writer.write('\n')
    writer.close()

if __name__ == "__main__":
    main(sys.argv)
