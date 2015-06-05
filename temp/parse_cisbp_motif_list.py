#/usr/bin/python
import numpy

dir_proj = "/Users/KANG/proj_db_infer_pipe/resources/"
fn_in = dir_proj + "cisbp_all_motifs.txt"

motifs = numpy.loadtxt(fn_in, dtype=str, delimiter="\t")

n_motifs_per_file = 200
n_files = numpy.ceil(len(motifs)/n_motifs_per_file)

for i in numpy.arange(n_files):
	ind_start = i*n_motifs_per_file
	ind_end = min((i+1)*n_motifs_per_file, len(motifs))
	temp_motifs = motifs[ind_start:ind_end]

	fn_out = dir_proj + "cisbp_all_motifs/cisbp_motifs_" + str(int(i+1)) + ".txt"
	numpy.savetxt(fn_out, temp_motifs, fmt="%s", newline="\n", delimiter="\t")
