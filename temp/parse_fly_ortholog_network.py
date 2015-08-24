#!/usr/bin/python
import sys
import numpy

species = sys.argv[1]
fn_gene_ortho = sys.argv[2]
fn_dmel_network = sys.argv[3]
fn_dmel_rids = sys.argv[4]
fn_dmel_gids = sys.argv[5]
fn_ortho_network = species + '.network.adjmtr'
fn_ortho_rids = species + '.rids'
fn_ortho_gids = species + '.gids'

dmel_network = numpy.loadtxt(fn_dmel_network)
dmel_rids = numpy.loadtxt(fn_dmel_rids, dtype=str)
dmel_gids = numpy.loadtxt(fn_dmel_gids, dtype=str)

index_rids = []
ortho_rids = []
index_gids = []
ortho_gids = []
lines = open(fn_gene_ortho, 'r').readlines()
for i in range(5, len(lines)):
	line = lines[i].split('\t')
	if (line[6].startswith(species)):
		if (line[0] in dmel_rids):
			index_rids.append(numpy.where(dmel_rids == line[0])[0][0])
			ortho_rids.append(line[5])
		if (line[0] in dmel_gids):
			index_gids.append(numpy.where(dmel_gids == line[0])[0][0])
			ortho_gids.append(line[5])

ortho_network = numpy.empty([len(index_rids), len(index_gids)])
for i in range(len(index_rids)):
	ortho_network[i,:] = dmel_network[index_rids[i], index_gids]

print ortho_network.shape

numpy.savetxt(fn_ortho_network, ortho_network, fmt='%s', delimiter='\t')
numpy.savetxt(fn_ortho_rids, ortho_rids, fmt='%s')
numpy.savetxt(fn_ortho_gids, ortho_gids, fmt='%s')
