#!/usr/bin/python
import os.path

fn_in = '/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_promoter_seq/ortho_scer+spar+smik+skud+sbay+scas+sklu.fasta'
dir_out = '/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_promoter_seq/ortho_scer+spar+smik+skud+sbay+scas+sklu/'

seqs = {}
lines = open(fn_in, 'r').readlines()
for i in range(len(lines)/2):
	seqs[lines[i*2].strip().strip('>')] = lines[i*2+1].strip()

for gid in sorted(seqs.keys()):
	gid = gid.strip('>')
	if len(gid.split('.')) == 1:
		fn_out = dir_out + gid
		writer = open(fn_out, 'w')
		writer.write('>%s\n%s\n' % (gid, seqs[gid]))
		writer.close()
	else:
		gid_split = gid.split('.')
		fn_out = dir_out + gid_split[0]
		with open(fn_out, 'a') as writer:
			writer.write('>%s.%s\n%s\n' % (gid_split[1], gid_split[0], seqs[gid]))
