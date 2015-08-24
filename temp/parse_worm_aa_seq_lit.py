#/usr/bin/python

fn_aa_fasta = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_aa_seq/c_elegans.PRJNA13758.WS240.protein.fa'
fn_aa_fasta_lit = '/home/mblab/ykang/proj_db_infer_pipe/resources/worm_aa_seq/c_elegans.PRJNA13758.WS240.protein.lit.fa'

lines = open(fn_aa_fasta, 'r').readlines()
writer = open(fn_aa_fasta_lit, 'w')

for i in range(len(lines)):
	if lines[i].startswith('>'):
		if i > 0:
			writer.write('%s\t%s\t%s\t\n%s\n' % (tf[0], tf[1], tf[2], seq.strip('*')))
		tf = lines[i].split()[:3]
		seq = ''
	else:
		seq += lines[i].strip()

writer.close()

lines = open(fn_aa_fasta_lit, 'r').readlines()

