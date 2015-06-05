#!/usr/bin/python
fn_hitdata = 'resources/fly_aa_seq/dmel-all-translation-r6.04.filtered.fasta_hitdata.txt'
fn_bed = 'resources/fly_aa_seq/dmel-all-translation-r6.04.filtered.fasta.bed'

writer = open(fn_bed, 'w')
lines = open(fn_hitdata, 'r').readlines()

prev_q_acc = ''
prev_q_from = ''
prev_q_to = ''
dup_counter = 0

for i in range(8, len(lines)):
	line = lines[i].split()
	if line[3] == 'superfamily' or line[3] == 'multi-dom':
		query = line[2].strip('>')
		q_from = line[5]
		q_to = line[6]
		q_acc = line[9]
		# print i, prev_q_acc
		if q_acc == prev_q_acc and q_from != prev_q_from and q_to != prev_q_to:
			dup_counter += 1
			writer.write('%s\t%s\t%s\t%s\n' % (query, q_from, q_to, query+'.'+q_acc+'-'+str(dup_counter)))
		else:
			dup_counter = 0
			writer.write('%s\t%s\t%s\t%s\n' % (query, q_from, q_to, query+'.'+q_acc))
		prev_q_acc = q_acc
writer.close()

