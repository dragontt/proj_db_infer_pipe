#!/usr/bin/python
import numpy

dir_res = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_baranski_0.15m_fatbody/'
fn_rids = dir_res + 'rids.fb'
fn_gids = dir_res + 'gids.fb'
fn_expr = dir_res + 'data.singles.expr'
fn_regulators_expr = dir_res + 'data.singles.regulators.expr'

rids = numpy.loadtxt(fn_rids, dtype=str)
gids = numpy.loadtxt(fn_gids, dtype=str)

writer = open(fn_regulators_expr, 'w')
lines = open(fn_expr, 'r').readlines()
for i in range(len(rids)):
	index = numpy.where(gids == rids[i])[0][0]
	writer.write('%s' % lines[index])

writer.close()


