#! /bin/bash
# Compute sequence percent identity of one query gene's sequence to those of all ortholog genes. 

DIR_IN=$1
DIR_RIDS_LIST=$2
DIR_OUT=$3

for RIDS_FILE in $DIR_RIDS_LIST/*
do
	qsub -P long -l h_vmem=4G /home/mblab/ykang/proj_db_infer_pipe/temp/compt_ortho_seq_pid_function.sh $DIR_IN $RIDS_FILE $DIR_OUT
done 

# function compute_pid() {
# 	while read rid
# 	do
# 		FILE=$1/$rid
# 		FN_OUT=$3/$(basename $FILE)
# 		clustalo --infile $FILE --outfile ${FN_OUT}.out --seqtype dna --distmat-out ${FN_OUT}.pim --full --percent-id --force
# 	done < $2
# }
