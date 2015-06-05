#!/usr/bin/python
import numpy
import glob
from os.path import basename

dir_res = '/home/mblab/ykang/proj_db_infer_pipe/resources/'
dir_expr = dir_res + 'fly_expression/'
fn_id_table = '/home/mblab/ykang/proj_db_infer_pipe/temp/fly_cg_fbgn_table.txt'

fn_expr_out = dir_res + 'fly_network_zeke_singles_net_full/data.singles.expr'
fn_cond_out = dir_res + 'fly_network_zeke_singles_net_full/conditions.singles.environmental'
fn_cond_fbgn_out = dir_res + 'fly_network_zeke_singles_net_full/conditions.singles.fbgn'

fns = sorted(glob.glob(dir_expr + '*.expr'))
cond = [None] * 158
expr = numpy.empty((14146, 158))
for i in range(158):
	cond[i] = basename(fns[i])
	expr[:,i] = numpy.loadtxt(fns[i])

numpy.savetxt(fn_expr_out, expr, fmt='%.1f', delimiter=' ', newline='\n')
		
writer = open(fn_cond_out, 'w')
for i in range(158):
	writer.write('%s\n' % cond[i])
writer.close()

ids = numpy.loadtxt(fn_id_table, dtype=str, skiprows=0)

writer = open(fn_cond_fbgn_out, 'w')
for i in range(len(cond)):
        line = cond[i].split('-')
        if line[0] == '00000':
                writer.write('\n')
        else:
                cgid = 'CG' + line[0]
                index = numpy.where(ids[:,0] == cgid)[0][0]
		# print index
                fbgnid = ids[index,1]
                writer.write('%s\n' % fbgnid)
writer.close()

