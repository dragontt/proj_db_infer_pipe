#!/usr/bin/python
import sys
import numpy 

# parsed_ortho = sys.argv[1]
parsed_ortho = ['Drosophila_simulans','Drosophila_sechellia','Drosophila_yakuba','Drosophila_erecta','Drosophila_ananassae']

fn_dmel_genes = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_baranski_singles_net_full/gids.fb'
dmel_genes = numpy.loadtxt(fn_dmel_genes, dtype=str)

ortho_list = {'Drosophila_ananassae':'Dana', 'Drosophila_erecta':'Dere', 'Drosophila_grimshawi':'Dgri', 'Drosophila_mojavensis':'Dmoj', 'Drosophila_persimilis':'Dper', 'Drosophila_pseudoobscura':'Dpse', 'Drosophila_sechellia':'Dsec', 'Drosophila_simulans':'Dsim', 'Drosophila_virilis':'Dvir', 'Drosophila_willistoni':'Dwil', 'Drosophila_yakuba':'Dyak'}
fn_ortho_conversion = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_base/gene_orthologs_fb_2015_03.tsv'
dir_promot_seqs = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_promoter_seq/'
# fn_dmel_ortho = dir_promot_seqs + 'rsat_dmel+'+ ortho_list[parsed_ortho] +'_upstream_-2000_0.filtered.fasta'
fn_dmel_ortho = dir_promot_seqs + 'rsat_dmel+'+ '+'.join([ortho_list[parsed_ortho[i]] for i in range(len(parsed_ortho))]) +'_upstream_-2000_0.filtered.fasta'

ortho_conversion = {}
for ortho in ortho_list.values():
	ortho_conversion[ortho] = {}
lines = open(fn_ortho_conversion, 'r').readlines()
for i in range(5, len(lines)):
	line = lines[i].split('\t')
	ortho = line[6].split('\\')[0]
	ortho_conversion[ortho][line[5]] = line[0]

promot_seqs = []

print 'Dmel'
fn_dmel = dir_promot_seqs + 'rsat_dmel_upstream_-2000_0.filtered.fasta'
lines = open(fn_dmel, 'r').readlines()
for i in range(len(lines)/2):
	gene = lines[i*2].strip().split('>')[1]
	seq = lines[i*2+1]
	if gene in dmel_genes:
		promot_seqs.append([gene, seq])

# for ortho in ortho_list.keys():
# for ortho in [parsed_ortho]:
for ortho in parsed_ortho:
	print ortho
	fn_ortho = dir_promot_seqs + 'rsat_'+ ortho +'_upstream_-2000_0.filtered.fasta'
	lines = open(fn_ortho, 'r').readlines()
	for i in range(len(lines)/2):
		gene = ortho_conversion[ortho_list[ortho]][lines[i*2].strip().split('>')[1]]
		seq = lines[i*2+1]
		if gene in dmel_genes:
			# promot_seqs.append([gene, seq])
			promot_seqs.append([lines[i*2].strip().split('>')[1], seq])

promot_seqs = numpy.array(promot_seqs, dtype=str)
promot_seqs = promot_seqs[promot_seqs[:,0].argsort()]

print fn_dmel_ortho
writer = open(fn_dmel_ortho, 'w')
for i in range(len(promot_seqs)):
	writer.write('>%s\n%s' % (promot_seqs[i,0], promot_seqs[i,1]))
writer.close()

