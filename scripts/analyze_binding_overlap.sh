# qsub -P long -l h_vmem=4G -N 'binding_overlap' ../../../../scripts/analyze_binding_overlap.sh ../np.adjmtr

NETWORK=$1
MAXRANK=$2
NUMBINS=$3
DIR_ANALYSIS=$4

REGS=/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_holstege/rids
GENES=/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_holstege/gids
CHIP_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_physical/chip_net.txt
PWM_NET=/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_physical/motif_net.txt

ADJLST=${DIR_ANALYSIS}/../${NETWORK##*/}.txt
NAME_ANALYSIS=${NETWORK##*/}; NAME_ANALYSIS=${NAME_ANALYSIS%.adjmtr}

if [ -z $NUMBINS ]; then NUMBINS=10; fi
MINRANK=$(bc <<< "scale=3; $MAXRANK/$NUMBINS")
if [ -z $DIR_ANALYSIS ]; then DIR_ANALYSIS=.; fi

sed -i 's/ /\t/g' ${NETWORK}

echo -e "REGULATOR\tTARGET\tCONFIDENCE" > ${ADJLST}
$HOME/proj_db_infer_pipe/scripts/adjmtr2interactions.rb -a ${NETWORK} -r ${REGS} -c ${GENES} >> ${ADJLST}

Rscript ~ykang/proj_db_infer_pipe/scripts/analyze_binding_overlap.r ${ADJLST} ${CHIP_NET} ${PWM_NET} ${DIR_ANALYSIS}/analysis.${NUMBINS}bins.top${MINRANK}to${MAXRANK}k.${NAME_ANALYSIS}.txt ${MAXRANK} ${NUMBINS}

rm ${ADJLST}
