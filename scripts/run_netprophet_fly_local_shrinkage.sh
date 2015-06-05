# In constructing fly network using local shrinkage, there are more than 500 variables and set use.Gram=FALSE in netprophet/CODE/global.lars.regulators.r
dir_res=/home/mblab/ykang/proj_db_infer_pipe/resources/fly_network_zeke_singles_net_full
netprophet -l -t ${dir_res}/data.singles.expr -r ${dir_res}/data.singles.regulators.expr -a ${dir_res}/allowed.adj -p ${dir_res}/data.singles.pert.adj -d ${dir_res}/rderank.cuffdiff.singles.adj -g ${dir_res}/gids.fb -f ${dir_res}/rids.fb  
