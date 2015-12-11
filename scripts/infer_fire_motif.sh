#! /bin/bash

fn_tf_list=$1		# a list of tf names
fn_fasta=$2 		# promoter sequence file
dir_binned_expr=$3 	# directory of binned expression files
fasta_len=$4		# length of promoter sequence (e.g. 600 for yeast, 2200 for fly)

cd $HOME/usr/FIRE-1.1a/
export FIREDIR=`pwd`

while read -a line
do
	tf=${line[0]}
	echo "__@__PROCESSING TF: $tf"
	perl fire.pl --expfiles=${dir_binned_expr}/$tf --exptype=discrete --fastafile_dna=${fn_fasta} -k=7 --seqlen_dna=${fasta_len} --nodups=1 --dorna=0 --dodnarna=0
done < ${fn_tf_list}
