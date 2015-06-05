#!/usr/bin/python

fn_gids = 'resources/fly_network_baranski_singles_net_full/gids.fb'
fn_fasta = 'resources/fly_promoter_seq/rsat_dmel_upstream_-1000_+200.lit.fasta'
fn_fasta_filt = 'resources/fly_promoter_seq/rsat_dmel_upstream_-1000_+200.filtered.fasta'

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
		gene = temp[1].strip()
		if (not gene in genes) and gene in gids:
			genes.append(gene)
			writer.write('>%s\n%s' % (gene, lines[i*2+1]))
writer.close()
