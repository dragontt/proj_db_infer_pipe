#!/bin/bash

tf_list=$1

DIR_PROJ=$HOME/proj_db_infer_pipe
DIR_RES=$DIR_PROJ/resources
# DIR_OUTPUT=$DIR_PROJ/output/fly_cisbp_dbd_pid
DIR_OUTPUT=$DIR_PROJ/output/fly_network_dbd_pid
TF_LIST=$DIR_RES/fly_aa_seq/${tf_list}

counter=0
while read line; do
	counter=$[$counter +1]
	if [ -e $DIR_OUTPUT/$line ]; then
		echo "Processing $line ... Existed ... $counter ... Done"
	else
		echo -n "Processing $line ... New ... $counter"
		# align against CIS-BP database motif sequences
		# bash $HOME/proj_db_infer_pipe/scripts/compt_dbd_pid.sh $line $DIR_PROJ $DIR_RES/fly_aa_seq/dmel-all-translation-r6.04.filtered.dbd.fasta $DIR_RES/cisbp_aa_seq/cisbp.dbd.fasta $DIR_OUTPUT 1
		
		# align against self TF DBD sequences
		bash $HOME/proj_db_infer_pipe/scripts/compt_dbd_pid.sh $line $DIR_PROJ $DIR_RES/fly_aa_seq/dmel-all-translation-r6.04.filtered.dbd.fasta $DIR_RES/fly_aa_seq/dmel-all-translation-r6.04.filtered.dbd.fasta $DIR_OUTPUT 0
		echo " ... Done"
	fi
done < $TF_LIST

