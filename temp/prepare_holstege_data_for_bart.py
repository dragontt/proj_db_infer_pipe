import numpy

dir_proj = '/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_raw_holstege/'

rids = numpy.loadtxt(dir_proj + 'rids', dtype=str)
gids = numpy.loadtxt(dir_proj + 'gids', dtype=str)

expr_genes = numpy.loadtxt(dir_proj +'holstege_logfc.txt', dtype=str, skiprows=1, usecols=[0])
expr_conds = numpy.loadtxt(dir_proj +'conditions.sysname', dtype=str)
expr_logfc = numpy.transpose(numpy.loadtxt(dir_proj +'holstege_logfc.txt', skiprows=1, usecols=(range(1,len(expr_conds)+1))))

gid_indices = []
for i in range(len(gids)):
	if gids[i] in expr_genes:
		gid_indices.append((i, numpy.where(expr_genes == gids[i])[0][0]))

# deal with replicates of conditions
unique_conds = numpy.unique(expr_conds)
for unique_cond in unique_conds:
	cond_idx = numpy.where(expr_conds==unique_cond)[0]
	if len(cond_idx) > 1:
		for i in range(len(cond_idx)):
			expr_conds[cond_idx[i]] = expr_conds[cond_idx[i]] +'_rep_'+ str(i) 

# data log fold change
expr_fc = numpy.ones((len(expr_conds), len(gids)))
for i in range(len(expr_conds)):
	for j in range(len(gid_indices)):
		gid_idx = gid_indices[j]
		expr_fc[i, gid_idx[0]] = numpy.exp2(expr_logfc[i, gid_idx[1]])
writer = open(dir_proj +'data.expr.fc.tsv', 'w')
for j in range(len(gids)):
	writer.write('\"%s\"\t' % gids[j])
writer.write('\n')
for i in range(len(expr_conds)):
	writer.write('\"%s"\t' % expr_conds[i])
	for j in range(len(gids)):
		writer.write('%.10f\t' % expr_fc[i,j])
	writer.write('\n')
writer.close()

# data perturbation
pert = numpy.transpose(numpy.loadtxt(dir_proj +'data.pert.adj'))
expr_pert = numpy.empty([len(expr_conds), len(gids)], dtype="S10")
for i in range(len(expr_conds)):
	for j in range(len(gids)):
		expr_pert[i,j] = 'TRUE' if pert[i,j] == 1 else 'FALSE'
writer = open(dir_proj +'data.pert.tsv', 'w')
for j in range(len(gids)):
	writer.write('\"%s\"\t' % gids[j])
writer.write('\n')
for i in range(len(expr_conds)):
	writer.write('\"%s"\t' % expr_conds[i])
	for j in range(len(gids)):
		writer.write('%s\t' % expr_pert[i,j])
	writer.write('\n')
writer.close()

# data regulator ids
writer = open(dir_proj +'rids.tsv', 'w')
for i in range(len(rids)):
	writer.write('\"%s\"\n' % rids[i])
writer.close()

