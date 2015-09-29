#!/usr/bin/python
import os.path

fn_ortho = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_base/gene_orthologs_fb_2015_03.tsv'
fn_in = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_promoter_seq/rsat_dmel+Dsim+Dsec+Dyak+Dere+Dana_upstream_-2000_+200.filtered.fasta'
dir_out = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_promoter_seq/rsat_dmel+Dsim+Dsec+Dyak+Dere+Dana_upstream_-2000_+200/'

seqs = {}
lines = open(fn_in, 'r').readlines()
for i in range(len(lines)/2):
	seqs[lines[i*2].strip()] = lines[i*2+1].strip()

lines = open(fn_ortho, 'r').readlines()
for i in range(5,len(lines)):
	line = lines[i].split()
	query_name = line[0]
	target_name = line[6].split('\\')[0] +'.'+ line[5]
	fn_out = dir_out + query_name
	if os.path.isfile(fn_out):
		with open(fn_out, 'a') as curr_file:
			target_seq_name = '>' + target_name.split('.')[1] 
			if target_seq_name in seqs.keys():
				curr_file.write('>%s\n%s\n' % (target_name, seqs[target_seq_name]))
	else:
		query_name = '>' + query_name
		if query_name in seqs.keys():
			curr_file = open(fn_out, 'w')
			curr_file.write('%s\n%s\n' % (query_name, seqs[query_name]))
			curr_file.close()
