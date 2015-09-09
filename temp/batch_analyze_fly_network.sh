#!/usr/bash
max_rank=$1
dir_analysis=$2

dir_output=/home/mblab/ykang/proj_db_infer_pipe/output
networks=()
networks+=(${dir_output}/fly_network_cellCycle_global_shrinkage/combined_model_full.adjmtr)
networks+=(${dir_output}/fly_network_cellCycle_motif_incorporated/cisbp_-2000_+200_fimo_known_motif/combined_network_np_motif_net_known_motif_resort.adjmtr)
networks+=(${dir_output}/fly_network_cellCycle_motif_incorporated/cisbp_-2000_+200_fimo_known_motif/combined_network_np_tf_merged_motif_net_known_motif_tf_merged_resort.adjmtr)
networks+=(${dir_output}/fly_network_cellCycle_motif_incorporated/fire_motifs_np_bin_20/combined_network_np_motif_net_fire_np_bin_20_resort.adjmtr)
networks+=(${dir_output}/fly_network_cellCycle_motif_incorporated/fire_motifs_np_tf_merged_dbd50_bin_20/combined_network_np_motif_net_fire_np_bin_20_resort_tf_merged.adjmtr)
networks+=(${dir_output}/fly_network_cellCycle_motif_incorporated/fire_ortho_motifs_np_bin_20/combined_network_np_ortho_motif_net_tf_merged_resort.adjmtr)
networks+=(${dir_output}/fly_network_cellCycle_motif_incorporated/fire_ortho_dmel+Dsim+Dsec_motifs_np_bin_20/combined_network_np_ortho_motif_net_tf_merged_resort.adjmtr)
networks+=(${dir_output}/fly_network_cellCycle_motif_incorporated/fire_ortho_dmel+Dsim+Dsec+Dyak+Dere+Dana_motifs_np_bin_20/combined_network_np_ortho_motif_net_tf_merged_resort.adjmtr)
networks+=(${dir_output}/fly_network_cellCycle_motif_incorporated/cisbp_-2000_+200_fimo_dbd_cutoff_cellCycle_np/combined_network_np_motif_net_dbd_cutoff_40.0_resort_0.adjmtr)
networks+=(${dir_output}/fly_network_cellCycle_motif_incorporated/cisbp_-2000_+200_fimo_dbd_cutoff_cellCycle_np_tf_merged_dbd50/combined_network_np_motif_net_dbd_cutoff_40_resort_tf_merged.adjmtr)

for network in ${networks[@]}
do
	# echo $(dirname $network)
	bash /home/mblab/ykang/proj_db_infer_pipe/scripts/analyze_compiled_chip_flynet_pwm.sh $network $max_rank $(dirname $network)/analysis_compiled_chip_flynet_pwm
done
