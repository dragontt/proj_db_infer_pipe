#! /usr/bin/python

import numpy
import matplotlib.pylab as plt

dir_proj = '/Users/KANG/cgscluster/proj_db_infer_pipe/'
fn_distmat = dir_proj + 'resources/crypto_aa_seq/clustal.dbd.distmat'
fn_distlst = dir_proj + 'resources/crypto_aa_seq/clustal.dbd.distlst'

dbds = numpy.loadtxt(fn_distmat, dtype=str, skiprows=1, usecols=[0])
dbds_adjmtr = numpy.loadtxt(fn_distmat, skiprows=1, usecols=range(1, len(dbds)+1))
print dbds_adjmtr.shape

# tfs_pair = []
# tfs_score = []
# for i in range(len(dbds)):
# 	print dbds[i]
# 	query = dbds[i].split(':')[0]
# 	for j in range(len(dbds)):
# 		target = dbds[j].split(':')[0]
# 		pair = query + '|' + target
# 		score = dbds_adjmtr[i,j]
# 		if not pair in tfs_pair:
# 			tfs_pair.append(pair)
# 			tfs_score.append(score)
# 		else:
# 			idx = tfs_pair.index(pair)
# 			tfs_score[idx] = score if score > tfs_score[idx] else tfs_score[idx]

# tfs_adjlst = zip(tfs_pair, tfs_score)
# tfs_adjlst_sorted = sorted(tfs_adjlst, key=lambda x: x[1])

# writer = open(fn_distlst, 'w')
# for pair in tfs_adjlst_sorted:
# 	pair_split = pair[0].split('|')
# 	if pair_split[0] != pair_split[1]:
# 		writer.write('%s\t%s\t%.2f\n' % (pair_split[0], pair_split[1], pair[1]))
# writer.close()

# find duplicated tfs indices
tfs_dup = []
for i in range(len(dbds)):
	tfs_dup.append(dbds[i].split(':')[0])
tfs_dup = numpy.array(tfs_dup)
tfs = numpy.unique(tfs_dup)
index_group = []
for i in range(len(tfs)):
	indices = numpy.where(tfs_dup == tfs[i])
	index_group.append(indices)

# sequeeze rows and columns
temp_adjmtr = numpy.max(dbds_adjmtr[index_group[0][0],:], axis=0)
for i in range(1, len(index_group)):
	max_row = numpy.max(dbds_adjmtr[index_group[i][0],:], axis=0)
	temp_adjmtr = numpy.vstack((temp_adjmtr, max_row))
temp_adjmtr = numpy.transpose(temp_adjmtr)
tfs_adjmtr = numpy.max(temp_adjmtr[index_group[0][0],:], axis=0)
for i in range(1, len(index_group)):
	max_row = numpy.max(temp_adjmtr[index_group[i][0],:], axis=0)
	tfs_adjmtr = numpy.vstack((tfs_adjmtr, max_row))

# # convert adjmtr to adjlst
# tfs_pair = []
# tfs_score = []
# for i in range(len(tfs)):
# 	for j in range(i+1, len(tfs)):
# 		pair = tfs[i] + ':' + tfs[j]
# 		score = tfs_adjmtr[i,j]
# 		tfs_pair.append(pair)
# 		tfs_score.append(score)

# tfs_adjlst = zip(tfs_pair, tfs_score)
# tfs_adjlst_sorted = sorted(tfs_adjlst, key=lambda x: x[1], reverse=True)

# writer = open(fn_distlst, 'w')
# for pair in tfs_adjlst_sorted:
# 	writer.write('%s\t%.2f\n' % (pair[0], pair[1]))
# writer.close()

# plot adjmtr
plt.imshow(tfs_adjmtr, interpolation='nearest', cmap=plt.cm.jet, extent=(0.5,10.5,0.5,10.5))
plt.colorbar()
plt.show()
