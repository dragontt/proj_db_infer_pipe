#!/bin/bash 
#for f in $HOME/proj_db_infer_pipe/resources/yeast_network_holstege/rids_list_10per/*;
#do
#	IN_TF_LIST=$f;
#	IN_TF_PWM=$HOME/proj_db_infer_pipe/output/yeast_network_holstege_orthologs/ortho_scer+smik+skud+sbay+scas+sklu/np_bart_tf_merged_orthos_motifs_pfm;
#	IN_PROMOTERS=$HOME/proj_db_infer_pipe/resources/yeast_promoter_seq/ortho_scer+smik+skud+sbay+scas+sklu.fasta;
#	IN_BACKGROUND=$HOME/proj_db_infer_pipe/resources/cisbp_all_species_bg_freq/Saccharomyces_cerevisiae.bmf;
#	OUT_FIMO=$HOME/proj_db_infer_pipe/output/yeast_network_holstege_orthologs/ortho_scer+smik+skud+sbay+scas+sklu/np_bart_tf_merged_orthos_promoters_fimo;
#
#	qsub -P long -l h_vmem=4G -N 'fimo_yeast' $HOME/proj_db_infer_pipe/scripts/compt_fimo_scan.sh $IN_TF_LIST $IN_TF_PWM $IN_PROMOTERS $IN_BACKGROUND $OUT_FIMO;
#done

#for f in $HOME/proj_db_infer_pipe/resources/fly_network_baranski_singles_net_full/rids_list_25per/*;
#do
#	IN_TF_LIST=$f;
#	IN_TF_PWM=$HOME/proj_db_infer_pipe/output/fly_network_cellCycle_orthologs/motif_inference_dmel+Dsim+Dsec+Dyak+Dere+Dana/combined_np_tf_merged_motifs_pfm;
#	IN_PROMOTERS=$HOME/proj_db_infer_pipe/resources/fly_promoter_seq/rsat_dmel+Dsim+Dsec+Dyak+Dere+Dana_upstream_-2000_+200.filtered.fasta;
#	IN_BACKGROUND=$HOME/proj_db_infer_pipe/resources/cisbp_all_species_bg_freq/Drosophila_melanogaster.bmf;
#	OUT_FIMO=$HOME/proj_db_infer_pipe/output/fly_network_cellCycle_orthologs/motif_inference_dmel+Dsim+Dsec+Dyak+Dere+Dana/combined_np_tf_merged_orthos_promoters_fimo;
#	
#	qsub -P long -l h_vmem=6G -N 'fimo_fly' $HOME/proj_db_infer_pipe/scripts/compt_fimo_scan.sh $IN_TF_LIST $IN_TF_PWM $IN_PROMOTERS $IN_BACKGROUND $OUT_FIMO;
#done

for f in $HOME/proj_db_infer_pipe/resources/yeast_network_holstege/rids_list_10per/*;
do
        IN_TF_LIST=$f;
        IN_TF_PWM=$HOME/proj_db_infer_pipe/output/yeast_network_raw_holstege_np_motif_incorporated/scertf_motif/scertf_pfm;
        IN_PROMOTERS=$HOME/proj_db_infer_pipe/resources/yeast_promoter_seq/s_cerevisiae.promoters.fasta;
        IN_BACKGROUND=$HOME/proj_db_infer_pipe/resources/cisbp_all_species_bg_freq/Saccharomyces_cerevisiae.bmf;
        OUT_FIMO=$HOME/proj_db_infer_pipe/output/yeast_network_raw_holstege_np_motif_incorporated/scertf_motif/scertf_fimo;

        qsub -P long -l h_vmem=4G -N 'fimo_yeast' $HOME/proj_db_infer_pipe/scripts/compt_fimo_scan.sh $IN_TF_LIST $IN_TF_PWM $IN_PROMOTERS $IN_BACKGROUND $OUT_FIMO;
done
