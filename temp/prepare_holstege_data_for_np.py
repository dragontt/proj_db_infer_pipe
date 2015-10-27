import numpy

dir_proj = '/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_raw_holstege/'

rids = numpy.loadtxt(dir_proj + 'rids', dtype=str)
gids = numpy.loadtxt(dir_proj + 'gids', dtype=str)

# rids/gids
dir_raw = dir_proj + 'netprophet_data/'
regulators = numpy.loadtxt(dir_raw + 'regulators', dtype=str)
genes = numpy.loadtxt(dir_raw + 'genes', dtype=str)
rid_indices = []
for i in range(len(rids)):
	if rids[i] in regulators:
		rid_indices.append((i, numpy.where(regulators == rids[i])[0][0]))
gid_indices = []
for i in range(len(gids)):
	if gids[i] in genes:
		gid_indices.append((i, numpy.where(genes == gids[i])[0][0]))

# allowed.adj
allowed = numpy.ones((len(rids), len(gids)))
for rid_idx in range(len(rids)):
	gid_idx = numpy.where(gids == rids[rid_idx])[0][0]
	allowed[rid_idx, gid_idx] = 0
numpy.savetxt('allowed.adj', allowed, fmt='%d', delimiter='\t')

# conditions
conditions = numpy.loadtxt('conditions.sysname', dtype=str)

# data.pert.adj
pert = numpy.zeros((len(gids), len(conditions)))
for cond_idx in range(len(conditions)):
	if conditions[cond_idx] in gids:
		gid_idx = numpy.where(gids == conditions[cond_idx])[0][0]
		pert[gid_idx, cond_idx] = 1
numpy.savetxt('data.pert.adj', pert, fmt='%d', delimiter='\t')

# data.expr
data_expr_raw = numpy.loadtxt(dir_raw +'data.expr')
data_expr = numpy.zeros((len(gids), len(conditions)))
for i in range(len(gid_indices)):
	gid_idx = gid_indices[i]
	data_expr[gid_idx[0],] = data_expr_raw[gid_idx[1],]
numpy.savetxt('data.expr', data_expr, fmt='%.10f', delimiter='\t')

# rdata.expr
rdata_expr = numpy.zeros((len(gids), len(conditions)))
for i in range(len(rids)):
	gid_idx = numpy.where(gids == rids[i])[0][0]
	rdata_expr[i,] = data_expr[gid_idx,]
numpy.savetxt('rdata.expr', rdata_expr, fmt='%.10f', delimiter='\t')

# prior.adj
prior_raw = numpy.loadtxt(dir_raw +'prior.adj')
prior = numpy.zeros((len(rids), len(gids)))
for i in range(len(rid_indices)):
	for j in range(len(gid_indices)):
		rid_idx, gid_idx = rid_indices[i], gid_indices[j]
		prior[rid_idx[0], gid_idx[0]] = prior_raw[rid_idx[1], gid_idx[1]]
numpy.savetxt('prior.adj', prior, fmt='%.10f', delimiter='\t')

# r_designs.adj
designs_raw = numpy.loadtxt(dir_raw +'r_designs.adj')
designs = numpy.zeros((len(rids), len(gids)))
for i in range(len(rid_indices)):
	for j in range(len(gid_indices)):
		rid_idx, gid_idx = rid_indices[i], gid_indices[j]
		designs[rid_idx[0], gid_idx[0]] = designs_raw[rid_idx[1], gid_idx[1]]
numpy.savetxt('r_designs.adj', designs, fmt='%.10f', delimiter='\t')
