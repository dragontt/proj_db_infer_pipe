# /usr/bin/python
import numpy
import matplotlib.pyplot as plt

dir_out = '/Users/KANG/cgscluster/proj_db_infer_pipe/'
datasets = ['baranski_singles_net_full', 'cellCycle', 'flyAtlas', 'geneDoseChange2L', 'genotypeByDiet', 'gravityResponse', 'lifeHistoryTraits', 'postFastingOlfactory', 'sexAntagonistic', 'toxicogenomicsLead']

conds = []
lassos = []
barts = []
des = [] 

for i in range(len(datasets)):
	fn_cond = dir_out + 'resources/fly_network_' + datasets[i] + '/conditions.count'
	conds.append(numpy.loadtxt(fn_cond, skiprows=1))	
	fn_lasso = dir_out + 'output/fly_network_' + datasets[i] + '_global_shrinkage/analysis_chip_support.txt'
	lassos.append(numpy.loadtxt(fn_lasso, skiprows=1))
	fn_bart = dir_out + 'output/fly_network_' + datasets[i] + '_bart/analysis_chip_support.txt'
	barts.append(numpy.loadtxt(fn_bart, skiprows=1))
	if i < 2:
		fn_de = dir_out + 'resources/fly_network_' + datasets[i] + '/analysis_chip_support.de.txt'
		des.append(numpy.loadtxt(fn_de, skiprows=1))

colors = ['r', 'g', 'b']
labels = ['de', 'lasso', 'bart']
cutoffs = ['2k', '4k', '6k', '8k', '10k', '12k', '14k', '16k', '18k', '20k']

ind = 9

de = []
lasso = []
bart = []
for i in range(len(datasets)):
	if i < 2:
		de.append(numpy.sum(des[i][0, :(ind+1)])/(ind+1))
	lasso.append(numpy.sum(lassos[i][0, :(ind+1)])/(ind+1))
	bart.append(numpy.sum(barts[i][0, :(ind+1)])/(ind+1))

plt.figure(num=None, figsize=(15,5), dpi=80)
plt.subplot(1,3,1)
y = []
for i in range(len(datasets)):
	y.append(conds[i][0])
plt.scatter(y[:2], de, c=colors[0], s=100, alpha=.5, label=labels[0])
plt.scatter(y, lasso, c=colors[1], s=100, alpha=.5, label=labels[1])
plt.scatter(y, bart, c=colors[2], s=100, alpha=.5, label=labels[2])
plt.xlabel('# samples')
plt.ylabel('ChIP support at ' + cutoffs[ind])
plt.legend(loc="upper left")

plt.subplot(1,3,2)
y = []
for i in range(len(datasets)):
	y.append(conds[i][1])
plt.scatter(y[:2], de, c=colors[0], s=100, alpha=.5)
plt.scatter(y, lasso, c=colors[1], s=100, alpha=.5)
plt.scatter(y, bart, c=colors[2], s=100, alpha=.5)
plt.xlabel('# experiments')
plt.ylabel('ChIP support at ' + cutoffs[ind])

plt.subplot(1,3,3)
y = []
for i in range(len(datasets)):
	y.append(conds[i][2])
plt.scatter(y[:2], de, c=colors[0], s=100, alpha=.5)
plt.scatter(y, lasso, c=colors[1], s=100, alpha=.5)
plt.scatter(y, bart, c=colors[2], s=100, alpha=.5)
plt.xlabel('# perturbations')
plt.ylabel('ChIP support at ' + cutoffs[ind])

plt.show()
