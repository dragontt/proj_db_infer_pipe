# qsub -P long -l h_vmem=4G -N 'analyze_network' ../../../../scripts/analyze_network.sh ../np.adjmtr 40

NETWORK=$1
MAXRANK=$2
NUMBINS=$3
DIR_ANALYSIS=$4

REGS=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_baranski_singles_net_full/rids.fb
GENES=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_baranski_singles_net_full/gids.fb
CHIP_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_chip/chip_net_genome_wide.adjmtr
PWM_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_net/resources/networks/physical/motif_net.adjmtr

fn=${NETWORK##*/}
fn=${fn%.adjmtr}

if [ -z $NUMBINS ]; then NUMBINS=10; fi

if [ -z $DIR_ANALYSIS ]; then DIR_ANALYSIS=.; fi

Rscript ~ykang/proj_db_infer_pipe/scripts/analyze_network.r ${NETWORK} ${REGS} ${GENES} ${CHIP_NET} ${DIR_ANALYSIS}/analysis_chip_support.${NUMBINS}bins.top$((MAXRANK / 10))to${MAXRANK}k.${fn}.txt ${MAXRANK} ${NUMBINS}
Rscript ~ykang/proj_db_infer_pipe/scripts/analyze_network.r ${NETWORK} ${REGS} ${GENES} ${PWM_NET} ${DIR_ANALYSIS}/analysis_pwm_support.${NUMBINS}bins.top$((MAXRANK / 10))to${MAXRANK}k.${fn}.txt ${MAXRANK} ${NUMBINS}
