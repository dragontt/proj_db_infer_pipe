# qsub -P long -l h_vmem=4G -N 'analyze_network' ../../../../scripts/analyze_network.sh ../np.adjmtr 40

NETWORK=$1
MAXRANK=$2
# REGS=$2
# GENES=$3
# CHIP_NET=$4
# PWM_NET=$5
REGS=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_baranski_singles_net_full/rids.fb
GENES=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_baranski_singles_net_full/gids.fb
CHIP_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_chip/chip_net_genome_wide.adjmtr
PWM_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_net/resources/networks/physical/motif_net.adjmtr

fn=${NETWORK##*/}
fn=${fn%.adjmtr}

Rscript ~ykang/proj_db_infer_pipe/scripts/analyze_network.r ${NETWORK} ${REGS} ${GENES} ${CHIP_NET} analysis_chip_support.top$((MAXRANK / 10))to${MAXRANK}k.${fn}.txt ${MAXRANK}

Rscript ~ykang/proj_db_infer_pipe/scripts/analyze_network.r ${NETWORK} ${REGS} ${GENES} ${PWM_NET} analysis_pwm_support.top$((MAXRANK / 10))to${MAXRANK}k.${fn}.txt ${MAXRANK}

