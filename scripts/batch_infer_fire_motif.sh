#!/bin/bash 
# Infer FIRE motifs 10 bins
for f in $HOME/proj_db_infer_pipe/resources/yeast_network_holstege/rids_list/*;
do
	TF_LIST=$f;
	PROMOTER=$HOME/proj_db_infer_pipe/resources/yeast_promoter_seq/s_cerevisiae.promoters.fasta;
	DIR_BINNED_EXPR=$HOME/proj_db_infer_pipe/output/yeast_network_raw_holstege_np_motif_inference/fire_motifs_bin_10/np_bins;
	SEQ_LENGTH=600;

	qsub -P long -l h_vmem=8G -N 'fire_yeast_bin_10' $HOME/proj_db_infer_pipe/scripts/infer_fire_motif.sh $TF_LIST $PROMOTER $DIR_BINNED_EXPR $SEQ_LENGTH;
done

# Infer FIRE motifs 50 bins
for f in $HOME/proj_db_infer_pipe/resources/yeast_network_holstege/rids_list/*;
do
	TF_LIST=$f;
	PROMOTER=$HOME/proj_db_infer_pipe/resources/yeast_promoter_seq/s_cerevisiae.promoters.fasta;
	DIR_BINNED_EXPR=$HOME/proj_db_infer_pipe/output/yeast_network_raw_holstege_np_motif_inference/fire_motifs_bin_50/np_bins;
	SEQ_LENGTH=600;

	qsub -P long -l h_vmem=8G -N 'fire_yeast_bin_50' $HOME/proj_db_infer_pipe/scripts/infer_fire_motif.sh $TF_LIST $PROMOTER $DIR_BINNED_EXPR $SEQ_LENGTH;
done