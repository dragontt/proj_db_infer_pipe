#!/bin/bash

tf_list=$1

DIR_RES=$HOME/proj_db_infer_pipe/resources
DIR_OUT=$HOME/proj_db_infer_pipe/output/fly_cisbp_dbd_pid
list=$DIR_RES/fly_aa_seq/${tf_list}
# list=$DIR_RES/crypto_network_orig/regulator.gids
counter=0
while read line; do
	counter=$[$counter +1]
	if [ -e $DIR_OUT/$line ]; then
		echo "Processing $line ... Existed ... $counter ... Done"
	else
		echo -n "Processing $line ... New ... $counter"
		bash $HOME/proj_db_infer_pipe/scripts/compt_dbd_pid.sh $line $DIR_RES/fly_aa_seq/dmel-all-translation-r6.04.filtered.dbd.fasta $DIR_RES/cisbp_aa_seq/cisbp.dbd.fasta 
		echo " ... Done"
	fi
done < $list

# list=$HOME/proj_db_infer_pipe/resources/crypto_network_orig/regulator.gids
# counter=0
# while read line; do
# 	counter=$[$counter +1]
# 	echo -n "Processing $line ... $counter"
# 	bash compt_dbd_pid.sh $line ~/proj_motifcomparison/resources/yeast.dbd.fasta
# 	echo " ... Done"
# done < $list
