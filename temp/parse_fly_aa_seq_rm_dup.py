#!/usr/bin/python
fn_aa_seq = 'resources/fly_aa_seq/dmel-all-translation-r6.04.lit.fasta'
fn_aa_seq_nodup = 'resources/fly_aa_seq/dmel-all-translation-r6.04.nodup.fasta'
fn_aa_dup_list = 'resources/fly_aa_seq/dmel-all-translation-r6.04.dup_list.txt'

lines = open(fn_aa_seq, 'r').readlines()
writer_seq = open(fn_aa_seq_nodup, 'w')
writer_dup = open(fn_aa_dup_list, 'w')

prev_name = lines[0].strip().split('>')[1]
prev_seq = lines[1].strip()
writer_seq.write('>%s\n%s\n' % (prev_name, prev_seq))
writer_dup.write('%s\t' % prev_name)
for i in range(1,len(lines)/2):
	curr_name = lines[i*2].strip().split('>')[1]
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
