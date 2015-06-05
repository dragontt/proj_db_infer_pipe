network <- as.matrix(read.table('/home/mblab/ykang/proj_db_infer_pipe/output/fly_network_singles_net_full_spearmanr/network_spearmanr.adjmtr'))
cutoff <- abs(network[order(abs(network),decreasing=TRUE)][40000])
cat('Abs edge score cutoff: ', cutoff)
network[which(abs(network) < cutoff)] <- 0
write.table(network, file='/home/mblab/ykang/proj_db_infer_pipe/output/fly_network_singles_net_full_spearmanr/network_spearmanr_40k.adjmtr', col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')

