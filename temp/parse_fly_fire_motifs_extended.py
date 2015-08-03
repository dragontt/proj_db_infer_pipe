#!/usr/bin/python
import sys
import numpy 

fn_curr = sys.argv[1]
fn_prev = sys.argv[2]
fn_extd = sys.argv[3]

motifs_curr = numpy.loadtxt(fn_curr, dtype=str)
motifs_prev = numpy.loadtxt(fn_prev, dtype=str)

motifs_diff = numpy.setdiff1d(motifs_prev[:,0], motifs_curr[:,0])
motifs_inds = []
for motif_diff in motifs_diff:
	motifs_inds.append(numpy.where(motifs_prev == motif_diff)[0][0])

motifs_extd = numpy.vstack((motifs_curr, motifs_prev[motifs_inds,:]))
numpy.savetxt(fn_extd, motifs_extd, fmt='%s', delimiter='\t')
