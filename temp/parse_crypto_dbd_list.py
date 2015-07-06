#/usr/bin/python
import numpy

dir_proj = "/home/mblab/ykang/proj_db_infer_pipe/resources/crypto_aa_seq/"
fn_in = dir_proj + "crNeoH99.peptides.lit.regulator.dbd.fasta"

# pids = numpy.loadtxt(fn_in, dtype=str, delimiter="\t")
pids = []
lines = open(fn_in, 'r').readlines()
for i in range(len(lines)/2):
	pids.append(lines[i*2].strip().split('>')[1])

n_pids_per_file = 50
n_files = int(numpy.ceil(len(pids)/n_pids_per_file))

for i in numpy.arange(n_files):
	ind_start = i*n_pids_per_file
	ind_end = min((i+1)*n_pids_per_file, len(pids))
	temp_pids = pids[ind_start:ind_end]

	fn_out = dir_proj + "pids_dbd_list/" + str(int(i+1)) + ".txt"
	numpy.savetxt(fn_out, temp_pids, fmt="%s", newline="\n", delimiter="\t")
