#!/usr/bin/python
import sys
import os
import glob
import numpy

fn_gene_ortho = sys.argv[1]
dir_dmel_network = sys.argv[2]
dir_combined_network = sys.argv[3]
dir_ortho_network = sys.argv[4]
ortho_names = sys.argv[5:]

if not dir_dmel_network.endswith('/'):
	dir_dmel_network += '/'
if not dir_combined_network.endswith('/'):
	dir_combined_network += '/'
if not dir_ortho_network.endswith('/'):
	dir_ortho_network += '/'

if not os.path.isdir(dir_combined_network):
	os.system('mkdir '+ dir_combined_network)

for i in range(len(ortho_names)):
	species = ortho_names[i]
	dir_converted = dir_ortho_network + species +'_dmel_converted/'
	os.system('mkdir '+ dir_converted)

	dict_conversion = {}
	lines = open(fn_gene_ortho, 'r').readlines()
	for i in range(5, len(lines)):
		line = lines[i].split('\t')
		if (line[6].startswith(species)):
			dict_conversion[line[5]] = line[0]

	fns_network_scores = glob.glob(dir_ortho_network + species + '/*')
	for fn_network_scores in fns_network_scores:
		fn_tf = os.path.basename(fn_network_scores)
		os.system('tail -n +2 '+ fn_network_scores +' > '+ dir_converted + dict_conversion[fn_tf])

fns_dmel_network_scores = glob.glob(dir_dmel_network +'*')
for fn_dmel_network_scores in fns_dmel_network_scores:
	rid_dmel = os.path.basename(fn_dmel_network_scores)
	cmdln = 'cat '+ fn_dmel_network_scores
	for species in ortho_names:
		fn_ortho = dir_ortho_network + species +'_dmel_converted/'+ rid_dmel
		if os.path.isfile(fn_ortho):
			cmdln += ' '+ fn_ortho
	cmdln += ' > '+ dir_combined_network + rid_dmel
	os.system(cmdln)

for i in range(len(ortho_names)):
	os.system('rm -r '+ dir_ortho_network + ortho_names[i] +'_dmel_converted')

