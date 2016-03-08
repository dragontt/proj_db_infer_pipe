#!/usr/bin/python
"""
Mask UCSC phastCons data to binary conserved elements with threshold.
"""
import glob
import sys

threshold = float(sys.argv[1])
print "Threshold =", threshold

# phastCons file directories
#dir_phastcons = '/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_promoter_seq/ucsc_phastCons/'
#fns = sorted(glob.glob(dir_phastcons + 'sacCer3.*.wigFixed'))
dir_phastcons = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_promoter_seq/ucsc_phastCons/'
fns = sorted(glob.glob(dir_phastcons + 'dm6.27way.phastCons.wigFix'))

# list of conserved hits in tuples (chr, start, end)
conserved_loci = []
# parse each phastCons file
for fn in fns:
	print "Parsing", fn
	
	lines = open(fn, 'r').readlines()
	for i in range(len(lines)):

		# get chromosome and start position annotation
		if lines[i].startswith('fixedStep'):
			chrom = lines[i].split()[1].split('=')[1]
			start = int(lines[i].split()[2].split('=')[1])
			count = 0
			flag = False
		# parse score
		else:
			score = float(lines[i].strip())
			if (score >= threshold) and (not flag):
				hit = (chrom, start+count)
				flag = True
			elif (score <= threshold) and (flag):
				hit += (start+count-1,)
				flag = False 
				conserved_loci.append(hit)
			count += 1

# save as conserved element like data
fn_output = dir_phastcons + '../ucsc_phastCons_threshould_'+ str(threshold) +'.txt'
writer = open(fn_output, 'w')
for i in range(len(conserved_loci)):
	writer.write('\t%s\t%d\t%d\n' % conserved_loci[i])
writer.close()
