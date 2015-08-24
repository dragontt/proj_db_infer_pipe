# use sequence id
import numpy

fn_seq = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_aa_seq/c_elegans.PRJNA13758.WS240.protein.nodup.fa'
fn_seq_filterd = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_aa_seq/c_elegans.PRJNA13758.WS240.protein.filtered.fa'
fn_rids = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_network_raw/wormbook.ce.tf_list.txt'

rids = numpy.loadtxt(fn_rids, dtype=str, delimiter='\t', skiprows=1, usecols=[0])

lines = open(fn_seq, 'r').readlines()
writer = open(fn_seq_filterd, 'w')
for i in range(len(lines)/2):
	gid = lines[i*2].strip().strip('>')
	if (gid in rids):
		writer.write('>%s\n%s' % (gid, lines[i*2+1]))
writer.close()
