# qsub -P long -l h_vmem=4G -N 'binding_overlap' ../../../../scripts/analyze_binding_overlap.sh ../np.adjmtr

NETWORK=$1
REGS=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_baranski_singles_net_full/rids.fb
GENES=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_baranski_singles_net_full/gids.fb
CHIP_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_chip/chip_net_genome_wide.adjlst
PWM_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_pwm/motif_net.adjlst

sed -i 's/ /\t/g' ${NETWORK}

echo -e "REGULATOR\tTARGET\tCONFIDENCE" > ${NETWORK}.txt
$HOME/proj_db_infer_pipe/scripts/adjmtr2interactions.rb -a ${NETWORK} -r ${REGS} -c ${GENES} >> ${NETWORK}.txt

R --no-save --slave --args ${NETWORK}.txt ${CHIP_NET} ${PWM_NET} < ~ykang/proj_db_infer_pipe/scripts/analyze_binding_overlap.r

rm ${NETWORK}.txt
