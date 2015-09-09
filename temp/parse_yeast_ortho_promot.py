#!/usr/python
import glob

dir_promot = '/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_promoter_seq/'
fn_gids = '/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_holstege/gids'
fn_name_conversion = dir_promot + 'Six_Species.name_conversion.txt'
fn_scerevisiae_seq = dir_promot + 's_cerevisiae.promoters.fasta'
dir_ortho_seq = dir_promot + 'Six_Species_Promoters/'

ortho_species = ['smik', 'skud', 'sbay']
fn_output = dir_promot + 'ortho_scer+' + '+'.join(ortho_species) + '.seq' 

name_conversion = {}
data = open(fn_name_conversion, 'r').readlines()
for i in range(len(data)):
	line = data[i].strip().split()
	name_conversion[line[2]] = line[1]

gids = []
data = open(fn_gids, 'r').readlines()
for i in range(len(data)):
	gids.append(data[i].strip())

scerevisiae_seqs = {}
data = open(fn_scerevisiae_seq, 'r').readlines()
for i in range(len(data)/2):
	if data[i*2].startswith('>'):
		name = data[i*2].strip().strip('>')
		if name in gids:
			seq = data[i*2+1].strip()
			scerevisiae_seqs[name] = seq

ortho_seqs = []
fns_ortho_seq = glob.glob(dir_ortho_seq + '*.seq')
for fn_ortho_seq in fns_ortho_seq:
	# print fn_ortho_seq
	data = open(fn_ortho_seq, 'r').readlines()
	gid = data[0].split()[0].strip('>')
	
	if gid in gids:
		if gid in scerevisiae_seqs.keys():
			ortho_seqs.append('>' + gid + '\n' + scerevisiae_seqs[gid] + '\n')
		
		for i in range(1,len(data)):
			if data[i].startswith('>'):
				species_num = data[i].split()[0].strip('>').split('.')
				species_name = name_conversion[species_num[0]]

				if species_name in ortho_species:
					if len(species_num) > 1:
						temp_name = gid + '.' + species_name + '.' + species_num[1]
					else:
						temp_name = gid + '.' + species_name
					temp_seq = data[i+1].strip()[-600:]
					ortho_seqs.append('>' + temp_name + '\n' + temp_seq + '\n')

writer = open(fn_output, 'w')
for i in range(len(ortho_seqs)):
	writer.write('%s' % ortho_seqs[i])
writer.close()
