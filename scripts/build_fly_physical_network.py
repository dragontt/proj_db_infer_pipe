#!/usr/bin/python
import numpy
import os.path

def main():
	""" Build gold standard from flyNet physical networks """
	#dir_flynet = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_net/resources/networks/physical/'
	#fn_rids = dir_flynet + 'rids.fb'
	#fn_gids = dir_flynet + 'gids.fb'
	#fn_chip_adjlst = dir_flynet + 'chip_net.txt'
	#fn_pwm_adjlst = dir_flynet + 'motif_net.txt'
	#fn_chip_adjmtr = dir_flynet + 'chip_net.adjmtr'
	#fn_pwm_adjmtr = dir_flynet + 'motif_net.adjmtr'
	
	# print "Convert ChIP adjmtr"
	# map_adjlst2adjmtr(fn_chip_adjlst, fn_chip_adjmtr, fn_rids, fn_gids, 1)
	# print "Convert PWM adjmtr"
	# map_adjlst2adjmtr(fn_pwm_adjlst, fn_pwm_adjmtr, fn_rids, fn_gids, 1)

	""" Build new chip network """
	# dir_chip = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_chip/'
	# fn_rids = dir_chip + 'rids.fb'
	# fn_gids = dir_chip + 'gids.fb'
	# fn_chip_adjlst = dir_chip + 'binding_interactions_only_genome_wide'
	# fn_chip_adjmtr = dir_chip + 'chip_net_only_genome_wide.adjmtr'
	
	# map_adjlst2adjmtr(fn_chip_adjlst, fn_chip_adjmtr, fn_rids, fn_gids, 0)

	""" Build new pwm network from CISBP direct evidence motifs """
	dir_proj = '/home/mblab/ykang/proj_db_infer_pipe/'
	fn_tf2motif = dir_proj + 'resources/cisbp_1.01/cisbp_tfs_Drosophila_melanogaster.txt'
	dir_fimo = dir_proj + 'output/fly_cisbp_-2000_+200_fimo_full/'
	fn_adjlst = dir_proj + 'resources/fly_physical_network_pwm/motif_net.adjlst'

	compile_adjlst(fn_tf2motif, dir_fimo, fn_adjlst)

def map_adjlst2adjmtr(fn_adjlst, fn_adjmtr, fn_rids, fn_gids, skiprows):
	rids = numpy.loadtxt(fn_rids, dtype=str)
	gids = numpy.loadtxt(fn_gids, dtype=str)
	lines = open(fn_adjlst, 'r').readlines()
	adjmtr = numpy.zeros((len(rids), len(gids)), dtype=numpy.int)	
	for i in range(skiprows,len(lines)):
		line = lines[i].split()
		r_index = numpy.where(rids == line[0])[0]
		g_index = numpy.where(gids == line[1])[0]
		if len(r_index) != 0 and len(g_index) != 0:
			adjmtr[r_index[0], g_index[0]] = 1

	numpy.savetxt(fn_adjmtr, adjmtr, fmt='%d', delimiter=' ')	

def compile_adjlst(fn_tf2motif, dir_fimo, fn_adjlst):
	tf2motif = numpy.loadtxt(fn_tf2motif, dtype=str, skiprows=1)
	writer = open(fn_adjlst, 'w')
	writer.write('REGULATOR\tTARGET\tSUM\tSUMP\tMAX\tMAXP\tCOUNT\tCOUNTP\n')
	for i in range(len(tf2motif)):
		[tf, motif] = tf2motif[i]
		fn_motif = dir_fimo + motif + '.summary'
		if os.path.isfile(fn_motif):
			lines = open(fn_motif, 'r').readlines()
			for line in lines:
				temp = line.split()
				writer.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (tf,temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],temp[7]))
	writer.close()

def parse_id_index(lst, ids):
	out_dict = {}
	for name in lst:
		index = numpy.where(ids == name)[0]
		if len(index) != 0:
			out_dict[name] = index[0]
	return out_dict

if __name__ == "__main__":
	main()
