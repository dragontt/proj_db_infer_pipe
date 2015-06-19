#!/usr/bin/python
import sys
import numpy

name = sys.argv[1]

dir_res = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_'+ name +'/'
# fn_pert = dir_res + 'data.pert'
# fn_fc = dir_res + 'data.expr.fc.tsv'
# fn_ispert = dir_res + 'data.ispert.tsv'
# fn_rids = dir_res + 'rids.fb'
fn_pert = dir_res + 'data.singles.pert.adj'
fn_fc = dir_res + 'data.singles.expr.fc.tsv'
fn_ispert = dir_res + 'data.singles.ispert.tsv'
fn_rids = dir_res + 'rids.fb'

row_header = open(fn_fc,'r').readline()
col_header = numpy.loadtxt(fn_fc, dtype=str, usecols=[0], skiprows=1, delimiter='\t')
x = numpy.transpose(numpy.loadtxt(fn_pert))

writer = open(fn_ispert, 'w')
writer.write('%s' % row_header)
for i in range(x.shape[0]):
	writer.write('%s\t' % col_header[i])
	for j in range(x.shape[1]):
		if x[i,j] == 1:
			#print i, j
			writer.write('TRUE\t')
		else:
			writer.write('FALSE\t')
	writer.write('\n')
writer.close() 

rids = numpy.loadtxt(fn_rids, dtype=str)
writer = open(fn_rids+'.tsv', 'w')
for i in range(len(rids)):
	writer.write('\"%s\"\n' % rids[i])
writer.close()
