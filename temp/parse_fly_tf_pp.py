#!/usr/bin/python
import numpy

fn_rids = 'resources/fly_netowrk_zeke_singles_net_full/rids.fb'
fn_conv = 'resources/fly_base/fbgn_fbtr_fbpp_fb_2015_01.tsv'
fn_pid_dup = 'resources/fly_aa_seq/dmel-all-translation-r6.04.dup_list.txt'
fn_fasta = 'resources/fly_aa_seq/dmel-all-translation-r6.04.nodup.fasta'

tfs = numpy.loadtxt(fn_rids, dtype=str)

pp = []
lines = open(fn_conv, 'r').readlines()
for i in range(len(lines)):
	if not (lines[i].startswith('##') or lines[i].startswith('\n')):
		tmp = lines[i].split()
		if len(tmp) > 2 and tmp[0] in tfs:
			pp.append(tmp[2])
pp = set(pp)

lines = open(fn_pid_dup, 'r').readlines()
for line in lines:
	temp = line.strip().split()
	if len(temp) > 1:
		for j in range(1,len(temp)):
			if temp[j] in pp:
				pp.remove(temp[j]) 

writer = open('resources/fly_aa_seq/dmel-all-translation-r6.04.filtered.fasta', 'w')
lines = open(fn_fasta, 'r').readlines()
for i in range(len(lines)/2):
	name = lines[i*2].strip().split('>')[1]
	seq = lines[i*2+1]
	if name in pp:
		writer.write('>%s\n%s' % (name, seq))
writer.close()

pp = sorted(list(pp))

writer = open('resources/fly_aa_seq/pids', 'w')
for i in range(len(pp)):
	writer.write('%s\n' % pp[i])
writer.close()

