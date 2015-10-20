#!/usr/bash
max_rank=$1
num_bins=$2

dir_output=/home/mblab/ykang/proj_db_infer_pipe/output
networks=()
networks+=(${dir_output}/yeast_network_holstege/np_bart_top400k.adjmtr)
networks+=(${dir_output}/yeast_network_holstege/np_bart_top400k_tf_merged_dbd50.adjmtr)
networks+=(${dir_output}/yeast_network_holstege_motif_incorporated/scertf_known_motif/combined_np_bart_tf_merged_motif_net_tf_merged.adjmtr)
networks+=(${dir_output}/yeast_network_holstege_motif_incorporated/fire_motifs_np_bart_tf_merged_dbd50_bin_20/combined_network_np_bart_tf_merged_net_fire_bin_20_tf_merged_resort.adjmtr)
networks+=(${dir_output}/yeast_network_holstege_motif_incorporated/fire_ortho_scer+spar_motifs_bin_20/combined_np_bart_tf_merged_motif_net_tf_merged.adjmtr)
networks+=(${dir_output}/yeast_network_holstege_motif_incorporated/cisbp_dbd_cutoff_holstege_np_bart_tf_merged_dbd50/combined_network_np_bart_motif_net_cisbp_dbd_cutoff_40_tf_merged_resort.adjmtr)

for network in ${networks[@]}
do
	# echo $(dirname $network)
	qsub -P long -l h_vmem=4G /home/mblab/ykang/proj_db_infer_pipe/scripts/analyze_binding_overlap.sh $network $max_rank $num_bins $(dirname $network)/analysis_binding_overlap
	# bash /home/mblab/ykang/proj_db_infer_pipe/scripts/analyze_binding_overlap.sh $network $max_rank $num_bins $(dirname $network)/analysis_binding_overlap
done
