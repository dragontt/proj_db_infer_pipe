# bash run_netprophet*.sh
dir_res=/home/mblab/ykang/proj_db_infer_pipe/resources/yeast_network_raw_holstege/netprophet_data
dir_out=/home/mblab/ykang/proj_db_infer_pipe/output/yeast_network_raw_holstege_np_global_test/

# force lars and rmpi version
#R CMD INSTALL /home/mblab/ykang/usr/netprophet_0.1/CODE/lars_0.9-8.tar.gz
#R CMD INSTALL /home/mblab/ykang/usr/netprophet_0.1/CODE/Rmpi_0.5-9.tar.gz

# RNA-seq
# netprophet -c -t ${dir_res}/data.singles.expr -r ${dir_res}/data.singles.regulators.expr -a ${dir_res}/allowed.adj -p ${dir_res}/data.singles.pert.adj -d ${dir_res}/rderank.cuffdiff.singles.adj -g ${dir_res}/gids.fb -f ${dir_res}/rids.fb 

# microarray
netprophet -m -c -t ${dir_res}/data.expr -r ${dir_res}/rdata.expr -a ${dir_res}/allowed.adj -p ${dir_res}/data.pert.adj -d ${dir_res}/prior.adj -g ${dir_res}/genes -f ${dir_res}/regulators -o ${dir_out} 
