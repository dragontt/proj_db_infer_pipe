fn_fasta = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_promoter_seq/rsat_celegans_upstream_-2000_+200.fasta'
fn_lit_fasta = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_promoter_seq/rsat_celegans_upstream_-2000_+200.lit.fasta'

lines = open(fn_fasta, 'r').readlines()
writer = open(fn_lit_fasta, 'w')
for i in range(len(lines)):
	if lines[i].startswith('>'):
		if i > 0:
			writer.write('>%s\n%s\n' % (gid, seq))
		gid = lines[i].split('|')[1]
		seq = ''
	else:
		if not lines[i].startswith('<'):
			seq += lines[i].strip()
writer.close()
