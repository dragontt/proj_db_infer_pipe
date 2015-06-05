# bash run_netprophet*.sh
DATASET=$1
dir_res=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_${DATASET}
dir_out=/home/mblab/ykang/proj_db_infer_pipe/output/fly_network_${DATASET}_global_shrinkage
# netprophet -c -t ${dir_res}/data.singles.expr -r ${dir_res}/data.singles.regulators.expr -a ${dir_res}/allowed.adj -p ${dir_res}/data.singles.pert.adj -d ${dir_res}/rderank.cuffdiff.singles.adj -g ${dir_res}/gids.fb -f ${dir_res}/rids.fb 
netprophet -m -c -t ${dir_res}/data.expr -r ${dir_res}/data.regulators.expr -a ${dir_res}/allowed.adjmtr -p ${dir_res}/data.pert -d ${dir_res}/rderank.cuffdiff.adj -g ${dir_res}/gids.fb -f ${dir_res}/rids.fb 
