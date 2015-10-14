#/usr/bin/python
fn_fasta = 'resources/fly_promoter_seq/rsat_dmel_upstream_-2000_0.fasta'
fn_fasta_lit = 'resources/fly_promoter_seq/rsat_dmel_upstream_-2000_0.lit.fasta'
fn_annot = 'resources/fly_base/fbgn_annotation_ID_fb_2015_01.tsv'

# Traverse the annotation file 3 times to get a full ID dictionary
annot1 = {}
annot2 = {}
annot3 = {}
lines = open(fn_annot, 'r').readlines()
for i in range(5, len(lines)):
	line = lines[i].split()
	if len(line) > 1:
		if len(line[0].split('\\')) > 1:
			line_1st = line[0].split('\\')[1]
		else:
			line_1st = line[0]
		annot1[line_1st.lower()] = line[1]
	if len(line) > 2:
		if not line[2].startswith('FBgn'):
			line_2nd = line[2]
		else:
			line_2nd = line[3]
		annot2[line_2nd.lower()] = line[1]
	if len(line) > 4:
		line_3rd = line[4].split(',')
		for j in range(len(line_3rd)):
			annot3[line_3rd[j].lower()] = line[1]

print 'Parsing annotation ... DONE'

# Get both gene symbol/ID and FBgn, and concatenate promoter sequence of each gene to one line
lines = open(fn_fasta, 'r').readlines()
writer = open(fn_fasta_lit, 'w')
missed = []
for i in range(len(lines)):
	if lines[i].startswith('>'):
		if i > 0:
			if tf.lower() in annot1.keys():
				writer.write('>%s|%s\n%s\n' % (tf, annot1[tf.lower()], seq.strip('*')))
			elif tf.lower() in annot2.keys():
				writer.write('>%s|%s\n%s\n' % (tf, annot2[tf.lower()], seq.strip('*')))
			elif tf.lower() in annot3.keys():
				writer.write('>%s|%s\n%s\n' % (tf, annot3[tf.lower()], seq.strip('*')))
			else:
				missed.append(tf)
				writer.write('>%s|\n%s\n' % (tf, seq.strip('*')))
		tf = lines[i].split()[0].split('|')[1].strip()
		seq = ''
	else:
		if not lines[i].startswith('<'):
			seq += lines[i].strip()
writer.close()
missed = list(set(missed))

print 'Parsing sequence ... DONE'
print 'Gene annotation missed: count =', len(missed)
print missed
