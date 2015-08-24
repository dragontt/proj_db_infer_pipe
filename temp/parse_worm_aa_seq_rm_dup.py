#!/usr/bin/python
fn_aa_seq = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_aa_seq/c_elegans.PRJNA13758.WS240.protein.lit.fa'
fn_aa_seq_nodup = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_aa_seq/c_elegans.PRJNA13758.WS240.protein.nodup.fa'
fn_aa_dup_list = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_aa_seq/c_elegans.PRJNA13758.WS240.protein.dup_list.txt'

lines = open(fn_aa_seq, 'r').readlines()
writer_seq = open(fn_aa_seq_nodup, 'w')
writer_dup = open(fn_aa_dup_list, 'w')

prev_name = lines[0].strip().strip('>').split('\t')[2]
prev_seq = lines[1].strip()
writer_seq.write('>%s\n%s\n' % (prev_name, prev_seq))
writer_dup.write('%s\t' % prev_name)
for i in range(1,len(lines)/2):
	curr_name = lines[i*2].strip().strip('>').split('\t')[2]
	curr_seq = lines[i*2+1].strip()
	if curr_seq == prev_seq:
		writer_dup.write('%s\t' % curr_name)
	else:
		writer_dup.write('\n%s\t' % curr_name)
		writer_seq.write('>%s\n%s\n' % (curr_name, curr_seq))
	prev_name = curr_name
	prev_seq = curr_seq

writer_seq.close()
writer_dup.close()
