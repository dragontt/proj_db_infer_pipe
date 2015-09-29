#!/usr/bin/python
import sys
import glob
import matplotlib.pyplot as plt

dir_in = sys.argv[1]
species = sys.argv[2]

orthos = {'yeast_index': ['4906', '6553', '4911', '476', '479'], 'yeast_name': ['smik','skud','sbay','scas','sklu'], 'fly': ['Dsim','Dsec','Dyak','Dere','Dana']}

# gids = numpy.loadtxt('/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_holstege/gids', dtype=str)
# gids = numpy.loadtxt('/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_baranski_singles_net_full/gids.fb', dtype=str)

pids = {}
for key in orthos[species]:
	pids[key] = []

fns_in = glob.glob(dir_in + '*.pim')
for fn_in in fns_in:
	lines = open(fn_in, 'r').readlines()
	for i in range(2,len(lines)):
		line = lines[i].split()
		pids[line[0].split('.')[0]].append(float(line[1]))

for key in pids.keys():
	l = pids[key]
	print key, reduce(lambda x, y: x + y, l) / len(l)

# yeast results:
# 4906 smik 50.067337
# 6553 skud 45.0745701945
# 4911 sbay 40.972110617
# 476 scas 25.5733251582
# 479 sklu 24.8165558631
#
# fly results:
# Dsim 64.6727991015
# Dsec 64.0593705472
# Dyak 52.2960690134
# Dere 53.8650128957
# Dana 31.4046373582
