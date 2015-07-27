#/usr/bin/python
import sys

species = sys.argv[1]
fn_fasta = sys.argv[2]
fn_fasta_lit = sys.argv[3]
fn_fbpp2fbgn = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_base/fbgn_fbtr_fbpp_fb_2015_01.tsv'
fn_ortho = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_base/gene_orthologs_fb_2015_03.tsv'

# Traverse the id conversion files to get a full dictionary
fbpp2fbgn = {}
lines = open(fn_fbpp2fbgn, 'r').readlines()
for i in range(6, len(lines)):
	line = lines[i].strip().split('\t')
	if len(line) == 3:
		fbpp2fbgn[line[2]] = line[0]

ortho = set()
lines = open(fn_ortho, 'r').readlines()
for i in range(5, len(lines)):
	if lines[i].split('\t')[6].startswith(species):
		gene = lines[i].split('\t')[5]
		ortho.add(gene)

print 'Parsing annotation ... DONE'

# Get both gene FBgn and FBpp, and concatenate promoter sequence of each gene to one line
lines = open(fn_fasta, 'r').readlines()
writer = open(fn_fasta_lit, 'w')
missed = set()
gene_count = 0
for i in range(len(lines)):
	if lines[i].startswith('>'):
		if i > 0:
			if (gene in fbpp2fbgn.keys()) and (fbpp2fbgn[gene] in ortho):
				writer.write('>%s|%s\n%s\n' % (fbpp2fbgn[gene], gene, seq.strip('*')))
				gene_count += 1
				if gene_count % 100 == 0:
					print gene_count
			else:
				missed.add(gene)
				# writer.write('>%s|\n%s\n' % (gene, seq.strip('*')))
		gene = lines[i].split('>')[1].split('|')[0]
		seq = ''
		# print gene, (gene in fbpp2fbgn.keys()), (fbpp2fbgn[gene] in ortho)
	else:
		if not lines[i].startswith('<'):
			seq += lines[i].strip()

writer.close()
missed = list(set(missed))

print 'Parsing sequence ... DONE'
print 'Gene annotation missed: count =', len(missed)
print missed
