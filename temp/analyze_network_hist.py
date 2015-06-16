#!/usr/bin/python
import numpy 
import pylab

dir_network = '/Users/KANG/cgscluster/proj_db_infer_pipe/output/'
# method = 'global_shrinkage'
# network = 'lasso.adjmtr'
method = 'bart'
network = 'fly_bart.adjmtr'

names = ['baranski_singles_net_full', 'cellCycle', 'flyAtlas', 'geneDoseChange2L', 'genotypeByDiet', 'gravityResponse', 'lifeHistoryTraits', 'postFastingOlfactory', 'sexAntagonistic', 'toxicogenomicsLead']

pylab.figure()

for i in range(len(names)):
	fn = dir_network + 'fly_network_' + names[i] + '_' + method + '/' + network
	x = numpy.loadtxt(fn)
	x = numpy.absolute(x[numpy.nonzero(x)])
	pylab.subplot(2,5,i+1)
	pylab.hist(x, 100, range=[0, 1.2], normed=1, histtype='bar')
	
pylab.show()
	
