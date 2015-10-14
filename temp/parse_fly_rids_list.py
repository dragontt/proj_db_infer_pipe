#/usr/bin/python
import sys
import os
import numpy

fn_in = sys.argv[1]
dir_out = sys.argv[2]

if not os.path.exists(dir_out):
	os.makedirs(dir_out)

rids = numpy.loadtxt(fn_in, dtype=str, delimiter="\t")

# n_rids_per_file = 25
n_rids_per_file = 1000
n_files = numpy.ceil(len(rids)/float(n_rids_per_file))

for i in numpy.arange(n_files):
	ind_start = i*n_rids_per_file
	ind_end = min((i+1)*n_rids_per_file, len(rids))
	temp_rids = rids[ind_start:ind_end]

	fn_out = dir_out + "/" + str(int(i+1)) + ".txt"
	numpy.savetxt(fn_out, temp_rids, fmt="%s", newline="\n", delimiter="\t")
