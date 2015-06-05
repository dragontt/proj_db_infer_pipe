#/usr/bin/python

fn_aa_fasta = 'resources/fly_base/dmel-all-translation-r6.04.fasta'
fn_aa_fasta_lit = 'resources/fly_aa_seq/dmel-all-translation-r6.04.lit.fasta'

lines = open(fn_aa_fasta, 'r').readlines()
writer = open(fn_aa_fasta_lit, 'w')

for i in range(len(lines)):
	if lines[i].startswith('>'):
		if i > 0:
			writer.write('>%s\n%s\n' % (tf, seq.strip('*')))
		tf = lines[i].split()[0].split('>')[1].strip()
		seq = ''
	else:
		seq += lines[i].strip()

writer.close()

lines = open(fn_aa_fasta_lit, 'r').readlines()

