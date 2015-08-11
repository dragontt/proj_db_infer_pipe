#!/usr/bin/python

ortho_list = ['Dana', 'Dere', 'Dgri', 'Dmoj', 'Dper', 'Dpse', 'Dsec', 'Dsim', 'Dvir', 'Dwil', 'Dyak']
fn_ortho_conversion = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_base/gene_orthologs_fb_2015_03.tsv'
dir_ortho = '/home/mblab/ykang/proj_db_infer_pipe/output/fly_network_cellCycle_orthologs/motif_inference/'

ortho_conversion = {}
for ortho in ortho_list:
	ortho_conversion[ortho] = {}

lines = open(fn_ortho_conversion, 'r').readlines()
for i in range(5, len(lines)):
	line = lines[i].split('\t')
	ortho = line[6].split('\\')[0]
	ortho_conversion[ortho][line[5]] = line[0]

dmel_motifs = {}
for ortho in ortho_list:
	lines = open(dir_ortho + ortho + '_motifs.txt', 'r').readlines()
	for i in range(len(lines)):
		line = lines[i].strip().split('\t')
		dmel_tf = ortho_conversion[ortho][line[0]]
		dmel_pwm = line[1]
		dmel_zscore = float(line[2])
		if not dmel_tf in dmel_motifs.keys():
			dmel_motifs[dmel_tf] = [dmel_pwm, dmel_zscore]
		else:
			if dmel_zscore > dmel_motifs[dmel_tf][1]:
				dmel_motifs[dmel_tf] = [dmel_pwm, dmel_zscore]

fn_output = dir_ortho + 'orthos_all_motifs.txt'
writer = open(fn_output, 'w')
for dmel_tf in sorted(dmel_motifs.keys()):
	writer.write('%s\t%s\t%s\n' % (dmel_tf, dmel_motifs[dmel_tf][0], str(dmel_motifs[dmel_tf][1])))
writer.close()


fn_dmel_np = '/home/mblab/ykang/proj_db_infer_pipe/output/fly_motif_inference/fire_7mers_cellCycle_np/motifs_np_bin_20.txt'
lines = open(fn_dmel_np, 'r').readlines()
for i in range(len(lines)):
	line = lines[i].strip().split('\t')
	dmel_tf = line[0]
	dmel_pwm = line[1]
	dmel_zscore= float(line[2])
	if not dmel_tf in dmel_motifs.keys():
		dmel_motifs[dmel_tf] = [dmel_pwm, dmel_zscore]
	else:
		if dmel_zscore >= dmel_motifs[dmel_tf][1]:
			dmel_motifs[dmel_tf] = [dmel_pwm, dmel_zscore]

fn_output = dir_ortho + 'combined_dmel_orthos_all_motifs.txt'
writer = open(fn_output, 'w')
for dmel_tf in sorted(dmel_motifs.keys()):
	writer.write('%s\t%s\t%s\n' % (dmel_tf, dmel_motifs[dmel_tf][0], str(dmel_motifs[dmel_tf][1])))
writer.close()
