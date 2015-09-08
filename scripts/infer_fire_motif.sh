#! /bin/bash

fn_fasta=$1 		# species' fasta file
dir_binned_expr=$2 	# directory of binned expression files
fn_tf_list=$3 		# a list of tf names
fasta_len=600
# fasta_len=2200

cd $HOME/usr/FIRE-1.1a/
export FIREDIR=`pwd`

while read -a line
do
	tf=${line[0]}
	echo "__@__PROCESSING TF: $tf"
	perl fire.pl --expfiles=${dir_binned_expr}/$tf --exptype=discrete --fastafile_dna=${fn_fasta} -k=7 --seqlen_dna=${fasta_len} --nodups=1 --dorna=0 --dodnarna=0
done < ${fn_tf_list}
