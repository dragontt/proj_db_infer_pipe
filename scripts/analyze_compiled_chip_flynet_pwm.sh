# qsub -P long -l h_vmem=4G -N 'binding_compiled' ../../../../scripts/analyze_compiled_chip_flynet_pwm.sh ../np.adjmtr 96.9 20 .

NETWORK=$1
MAXRANK=$2
NUMBINS=$3
DIR_ANALYSIS=$4

REGS=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_cellCycle/rids_full.fb
GENES=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_cellCycle/gids_full.fb

CHIP_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_chip/chip_net_genome_wide.adjmtr
PWM_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_pwm/flynet_motif_net.adjmtr

# CHIP_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_chip/chip_net_genome_wide.intersect_with_pwm.adjmtr
# PWM_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_pwm/flynet_motif_net.intersect_with_chip.adjmtr

fn=${NETWORK##*/}
fn=${fn%.adjmtr}

if [ -z $NUMBINS ]; then NUMBINS=10; fi
MINRANK=$(bc <<< "scale=3; $MAXRANK/$NUMBINS")
if [ -z $DIR_ANALYSIS ]; then DIR_ANALYSIS=.; fi

Rscript ~ykang/proj_db_infer_pipe/scripts/analyze_network.r ${NETWORK} ${REGS} ${GENES} ${CHIP_NET} ${DIR_ANALYSIS}/analysis_chip_support.${NUMBINS}bins.top${MINRANK}to${MAXRANK}k.${fn}.txt ${MAXRANK} ${NUMBINS}
Rscript ~ykang/proj_db_infer_pipe/scripts/analyze_network.r ${NETWORK} ${REGS} ${GENES} ${PWM_NET} ${DIR_ANALYSIS}/analysis_pwm_support.${NUMBINS}bins.top${MINRANK}to${MAXRANK}k.${fn}.txt ${MAXRANK} ${NUMBINS}

# Rscript ~ykang/proj_db_infer_pipe/scripts/analyze_network.r ${NETWORK} ${REGS} ${GENES} ${CHIP_NET} ${DIR_ANALYSIS}/analysis_chip_support.chip_pwm_intersected.${NUMBINS}bins.top${MINRANK}to${MAXRANK}k.${fn}.txt ${MAXRANK} ${NUMBINS}
# Rscript ~ykang/proj_db_infer_pipe/scripts/analyze_network.r ${NETWORK} ${REGS} ${GENES} ${PWM_NET} ${DIR_ANALYSIS}/analysis_pwm_support.chip_pwm_intersected.${NUMBINS}bins.top${MINRANK}to${MAXRANK}k.${fn}.txt ${MAXRANK} ${NUMBINS}
