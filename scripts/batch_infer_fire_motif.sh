#!/bin/bash 
# Infer scer FIRE motifs 20 bins
for f in $HOME/proj_db_infer_pipe/resources/yeast_network_holstege/rids_list/*;
do
	TF_LIST=$f;
	PROMOTER=$HOME/proj_db_infer_pipe/resources/yeast_promoter_seq/s_cerevisiae.promoters.fasta;
	DIR_BINNED_EXPR=$HOME/proj_db_infer_pipe/output/yeast_network_raw_holstege_np_bart_motif_incorporated/fire_scer/np_bart_bins/;
	SEQ_LENGTH=600;

	qsub -P long -l h_vmem=8G -N 'motif_yeast' $HOME/proj_db_infer_pipe/scripts/infer_fire_motif.sh $TF_LIST $PROMOTER $DIR_BINNED_EXPR $SEQ_LENGTH;
done

# Infer dmel FIRE motifs 20 bins
for f in $HOME/proj_db_infer_pipe/resources/fly_network_baranski_singles_net_full/rids_list_25per/*;
do
	TF_LIST=$f;
	PROMOTER=$HOME/proj_db_infer_pipe/resources/fly_promoter_seq/rsat_dmel_upstream_-2000_+200.filtered.fasta;
	DIR_BINNED_EXPR=$HOME/proj_db_infer_pipe/output/fly_network_cellCycle_np_bart_motif_incorporated/fire_dmel/np_bart_bins/;
	SEQ_LENGTH=2200;

	qsub -P long -l h_vmem=8G -N 'motif_fly' $HOME/proj_db_infer_pipe/scripts/infer_fire_motif.sh $TF_LIST $PROMOTER $DIR_BINNED_EXPR $SEQ_LENGTH;
done