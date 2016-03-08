#!/usr/bin/python
"""
Mask promoters using UCSC conserved elements.
"""

import sys
import os
import argparse
import numpy
import re

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Mask promoters using UCSC conserved elements.")
    parser.add_argument('-g', '-fn_gids', dest='fn_gids')
    parser.add_argument('-a', '-fn_gene_annot', dest='fn_gene_annot')
    parser.add_argument('-f', '-fn_gene_annot_format', dest='fn_gene_annot_format')
    parser.add_argument('-e', '-fn_conselem', dest='fn_conselem')
    parser.add_argument('-l', '-promoter_length', dest='promoter_length', type=int)
    parser.add_argument('-t', '-promoter_tss_shift', dest='promoter_tss_shift', type=int, default=0) # -2000_+200 fly promoter, promoter_tss_shift=200
    parser.add_argument('-o', '-fn_output_mask', dest='fn_output_mask')
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)

    # read files
    gids = numpy.loadtxt(parsed.fn_gids, dtype=str)
    cons_elem = numpy.loadtxt(parsed.fn_conselem, dtype=str, usecols=(1,2,3), delimiter='\t')
    prom_len = parsed.promoter_length
    tss_shift = parsed.promoter_tss_shift

    if parsed.fn_gene_annot_format == 'sgd':
        gene_annot = numpy.loadtxt(parsed.fn_gene_annot, dtype=str, usecols=(1,2,3,6,7), delimiter='\t')
    elif parsed.fn_gene_annot_format == 'flybase':
        gene_annot_raw = numpy.loadtxt(parsed.fn_gene_annot, dtype=str, usecols=(1,4), delimiter='\t', skiprows=6)
        gene_annot = []
        for i in range(len(gene_annot_raw)):
            if gene_annot_raw[i,1] != '':
                gene_id = gene_annot_raw[i,0]
                line = re.split('\:|\..|\(|\)', gene_annot_raw[i,1])
                chrm = 'chr' + line[0]
                sign = '+' if line[3] == '1' else '-'
                (start, end) = (line[1], line[2])
                gene_annot.append([gene_id, chrm, sign, start, end])
        gene_annot = numpy.array(gene_annot)
        print 'flybase dimension:',gene_annot.shape
    else: 
        print "No gene annotation format specified."
    
    # mask promoters 
    prom_mask = {}
    chrms = numpy.intersect1d( numpy.unique(gene_annot[:,1]), numpy.unique(cons_elem[:,0]) )
    print "chromosome used:", chrms
    for chrm in chrms:
        subset_gene_annot = gene_annot[numpy.where(gene_annot[:,1]==chrm)[0]][:,(0,2,3,4)]
        subset_cons_elem = cons_elem[numpy.where(cons_elem[:,0]==chrm)[0]][:,(1,2)]
        subset_cons_elem = numpy.array(subset_cons_elem, dtype=int)
        
        # mask promoter of each target gene
        for i in range(len(subset_gene_annot)):
            
            # direction of transcription is positive
            if subset_gene_annot[i,1] == '+': 
                (gid, tss) = (subset_gene_annot[i,0], int(subset_gene_annot[i,2])+tss_shift)
                if gid in gids:
                    (prom_start, prom_end) = (tss-prom_len, tss-1)
                    elem_lefts = subset_cons_elem[:,0]
                    elem_rights = subset_cons_elem[:,1]
                    # conserved regions within the promoter
                    inds = numpy.intersect1d(numpy.where(elem_lefts>=prom_start)[0],
                        numpy.where(elem_rights<=prom_end)[0])
                    prom_mask[gid] = subset_cons_elem[inds] - tss + prom_len
                    # conserved region intersect with promoter ends
                    index = numpy.intersect1d(numpy.where(elem_lefts<=prom_start)[0],
                        numpy.where(elem_rights>=prom_start)[0])
                    if index.size > 0:
                        prom_mask[gid] = numpy.vstack((prom_mask[gid], 
                            [0, subset_cons_elem[index[0],1] - tss + prom_len]))
                    index = numpy.intersect1d(numpy.where(elem_lefts<=prom_end)[0],
                        numpy.where(elem_rights>=prom_end)[0])
                    if index.size > 0:
                        prom_mask[gid] = numpy.vstack((prom_mask[gid], 
                            [subset_cons_elem[index[0],0] - tss + prom_len, prom_len]))

            # direction of transcription is negative
            elif subset_gene_annot[i,1] == '-': 
                (gid, tss) = (subset_gene_annot[i,0], int(subset_gene_annot[i,3])-tss_shift)
                if gid in gids:
                    (prom_start, prom_end) = (tss+prom_len, tss+1)
                    elem_lefts = subset_cons_elem[:,0]
                    elem_rights = subset_cons_elem[:,1]
                    # conserved regions within the promoter
                    inds = numpy.intersect1d(numpy.where(elem_rights<=prom_start)[0],
                        numpy.where(elem_lefts>=prom_end)[0])
                    prom_mask[gid] = numpy.fliplr(prom_len + tss - subset_cons_elem[inds])
                    # conserved region intersect with promoter ends
                    index = numpy.intersect1d(numpy.where(elem_rights>=prom_start)[0],
                        numpy.where(elem_lefts<=prom_start)[0])
                    if index.size > 0:
                        prom_mask[gid] = numpy.vstack((prom_mask[gid],
                            [0, prom_len + tss - subset_cons_elem[index[0],0]]))
                    index = numpy.intersect1d(numpy.where(elem_rights>=prom_end)[0],
                        numpy.where(elem_lefts<=prom_end)[0])
                    if index.size > 0:
                        prom_mask[gid] = numpy.vstack((prom_mask[gid],
                            [prom_len + tss - subset_cons_elem[index[0],1], prom_len]))

            else:
                print "Error in starin annotation."

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
