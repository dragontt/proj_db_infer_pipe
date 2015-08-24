# use sequence id
import numpy

fn_seq = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_promoter_seq/rsat_celegans_upstream_-2000_+200.lit.fasta'
fn_seq_filterd = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_promoter_seq/rsat_celegans_upstream_-2000_+200.filtered.fasta'
fn_gids = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_network_raw/wormbase.WS240.ce.gene_ids.txt'

gene_web_ids = numpy.loadtxt(fn_gids, dtype=str, delimiter='\t', skiprows=4, usecols=[2])
gene_pub_ids = numpy.loadtxt(fn_gids, dtype=str, delimiter='\t', skiprows=4, usecols=[3])
gene_seq_ids = numpy.loadtxt(fn_gids, dtype=str, delimiter='\t', skiprows=4, usecols=[4])

lines = open(fn_seq, 'r').readlines()
writer = open(fn_seq_filterd, 'w')
gids_parsed = set()
for i in range(len(lines)/2):
	gid = lines[i*2].strip().strip('>')
	if (not gid in gids_parsed) and (gid in gene_seq_ids):
		gid_parsed = gene_web_ids[numpy.where(gene_seq_ids == gid)[0][0]]
		writer.write('>%s\n%s' % (gid_parsed, lines[i*2+1]))
		gids_parsed.add(gid_parsed)
	elif (not gid in gids_parsed) and (gid in gene_pub_ids):
		gid_parsed = gene_web_ids[numpy.where(gene_pub_ids == gid)[0][0]]
		writer.write('>%s\n%s' % (gid_parsed, lines[i*2+1]))
		gids_parsed.add(gid_parsed)
writer.close()
