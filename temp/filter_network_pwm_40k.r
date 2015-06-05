network <- as.matrix(read.table('/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_pwm/motif_net.adjmtr'))
# cutoff <- abs(network[order(abs(network),decreasing=TRUE)][20000])
cutoff <- abs(network[order(abs(network),decreasing=TRUE)][100000])
cat('Abs edge score cutoff: ', cutoff)
network[which(abs(network) < cutoff)] <- 0
network[which(abs(network) > cutoff)] <- 1
# write.table(network, file='/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_pwm/motif_net_20k.adjmtr', col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
write.table(network, file='/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_pwm/motif_net_100k.adjmtr', col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
