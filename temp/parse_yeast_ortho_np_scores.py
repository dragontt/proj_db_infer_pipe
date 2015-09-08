#!/usr/python
import glob
import os

fn_gids = '/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_holstege/gids'
fn_seqs = '/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_promoter_seq/Six_Species_Promoters.seq'
dir_np_scores = '/home/mblab/ykang/proj_db_infer_pipe/output/yeast_network_holstege_orthologs/np_bart_tf_merged_bin_20/'
dir_combined_scores = '/home/mblab/ykang/proj_db_infer_pipe/output/yeast_network_holstege_orthologs/np_bart_tf_merged_orthos_bin_20/'

orthos = {}
data = open(fn_gids, 'r').readlines()
for i in range(len(data)):
	orthos[data[i].strip()] = []

data = open(fn_seqs, 'r').readlines()
for i in range(len(data)/2):
	name = data[i*2].strip('>').strip()
	orthos[name.split('.')[0]].append(name)

fns_np_scores = glob.glob(dir_np_scores + '*')
for fn_np_scores in fns_np_scores:
	tf = os.path.basename(fn_np_scores)
	scores = open(fn_np_scores, 'r').readlines()
	writer = open(dir_combined_scores + tf, 'w')
	for i in range(1,len(scores)):
		[name, score] = scores[i].strip().split()
		for j in range(len(orthos[name])):
			writer.write('%s\t%s\n' % (orthos[name][j], score))
	writer.close()
