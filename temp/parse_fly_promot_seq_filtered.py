#!/usr/bin/python
import sys

fn_gids = sys.argv[1]
fn_fasta = sys.argv[2]
fn_fasta_filt = sys.argv[3]

gids = []
lines = open(fn_gids, 'r').readlines()
for line in lines:
	gids.append(line.strip())

genes = []
lines = open(fn_fasta, 'r').readlines()
writer = open(fn_fasta_filt, 'w')
for i in range(len(lines)/2):
	temp = lines[i*2].split('>')[1].split('|')
	if len(temp) > 1:
		gene = temp[0].strip()
		# gene = temp[1].strip()
		if (not gene in genes) and (gene in gids):
			genes.append(gene)
			writer.write('>%s\n%s' % (gene, lines[i*2+1]))
writer.close()
