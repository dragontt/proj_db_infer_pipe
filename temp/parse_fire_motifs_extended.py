#!/usr/bin/python
import sys
import numpy 

fn_curr = sys.argv[1]
fn_prev = sys.argv[2]
fn_extd = sys.argv[3]

motifs_curr = numpy.loadtxt(fn_curr, dtype=str)
motifs_prev = numpy.loadtxt(fn_prev, dtype=str)

# # use current motifs, and append previously avail motifs as complement
# motifs_diff = numpy.setdiff1d(motifs_prev[:,0], motifs_curr[:,0])
# motifs_inds = []
# for motif_diff in motifs_diff:
# 	motifs_inds.append(numpy.where(motifs_prev == motif_diff)[0][0])

# motifs_extd = numpy.vstack((motifs_curr, motifs_prev[motifs_inds,:]))
# numpy.savetxt(fn_extd, motifs_extd, fmt='%s', delimiter='\t')

# take the union of prev. and curr. motifs based on mutual info z-score
motifs_union = {}
for i in range(len(motifs_prev)):
	motifs_union[motifs_prev[i,0]] = [motifs_prev[i,1], float(motifs_prev[i,2])]

counter_no_improvement = 0

for i in range(len(motifs_curr)):
	if motifs_curr[i,0] in motifs_union.keys():
		
		if float(motifs_curr[i,2]) > motifs_union[motifs_curr[i,0]][1]:
			motifs_union[motifs_curr[i,0]] = [motifs_curr[i,1], float(motifs_curr[i,2])]
		else:
			counter_no_improvement += 1
	
	else:
		motifs_union[motifs_curr[i,0]] = [motifs_curr[i,1], float(motifs_curr[i,2])]

if counter_no_improvement == len(motifs_curr):
	print "No motif improvement."
else:
	print "Motif improved: " + str(len(motifs_curr) - counter_no_improvement)
	writer = open(fn_extd, 'w')
	for motif in sorted(motifs_union.keys()):
		writer.write('%s\t%s\t%s\n' % (motif, motifs_union[motif][0], str(motifs_union[motif][1])))
	writer.close()
