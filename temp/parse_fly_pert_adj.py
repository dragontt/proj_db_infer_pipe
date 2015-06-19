#!/usr/bin/python 
import numpy

dir_res = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_baranski_0.15m_fatbody/'
fn_expr = dir_res + 'data.singles.expr'
fn_gids = dir_res + 'gids.fb'
fn_cond = dir_res + 'conditions.singles.fbgn'

expr = numpy.loadtxt(fn_expr)
gids = numpy.loadtxt(fn_gids, dtype=str)

pert = numpy.zeros(expr.shape)

lines = open(fn_cond, 'r').readlines()
for i in range(len(lines)):
	line = lines[i].strip()
	if line != '':
		index = numpy.where(gids == line)[0][0]
		# print index
		pert[index,i] = 1

numpy.savetxt(dir_res + 'data.singles.pert.adj', pert, fmt='%d', delimiter=' ', newline='\n')
	
