#!/usr/bin/python
import sys
import os

fn_inferred_motif = sys.argv[1]
dir_cisbp_pfm = sys.argv[2]
dir_inferred_pfm = sys.argv[3]

dir_cisbp_pfm += '/' if not dir_cisbp_pfm.endswith('/') else ''
dir_inferred_pfm += '/' if not dir_inferred_pfm.endswith('/') else ''

os.system('mkdir -p '+ dir_inferred_pfm)

lines = open(fn_inferred_motif, 'r').readlines()
for i in range(1, len(lines)):
	id_motif, id_cisbp = lines[i].split()[:2]
	os.system('ln -s ../'+ dir_cisbp_pfm + id_cisbp +' '+ dir_inferred_pfm + id_motif)
