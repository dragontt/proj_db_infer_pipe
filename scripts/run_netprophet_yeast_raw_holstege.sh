# bash run_netprophet*.sh
dir_res=/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_raw_holstege
dir_out=/home/mblab/ykang/proj_db_infer_pipe/output/fly_network_raw_holstege

# RNA-seq
# netprophet -c -t ${dir_res}/data.singles.expr -r ${dir_res}/data.singles.regulators.expr -a ${dir_res}/allowed.adj -p ${dir_res}/data.singles.pert.adj -d ${dir_res}/rderank.cuffdiff.singles.adj -g ${dir_res}/gids.fb -f ${dir_res}/rids.fb 

# microarray
netprophet -m -c -t ${dir_res}/data.expr -r ${dir_res}/rdata.expr -a ${dir_res}/allowed.adj -p ${dir_res}/data.pert.adj -d ${dir_res}/prior.adj -g ${dir_res}/gids -f ${dir_res}/rids 
