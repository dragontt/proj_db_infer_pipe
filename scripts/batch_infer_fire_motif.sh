#!/bin/bash 
for f in $HOME/proj_db_infer_pipe/resources/yeast_network_holstege/rids_list/*;
do
	TF_LIST=$f;
	PROMOTER=$HOME/proj_db_infer_pipe/resources/yeast_promoter_seq/s_cerevisiae.promoters.fasta;
	DIR_BINNED_EXPR=$HOME/proj_db_infer_pipe/output/yeast_motif_inference/fire_7mers_holstege_np_bart/np_bart_bin_20;
	SEQ_LENGTH=600;

	qsub -P long -l h_vmem=8G -N 'fire_yeast' $HOME/proj_db_infer_pipe/scripts/infer_fire_motif.sh $TF_LIST $PROMOTER $DIR_BINNED_EXPR $SEQ_LENGTH;
done
