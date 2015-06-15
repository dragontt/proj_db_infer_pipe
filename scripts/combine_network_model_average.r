# Combine NetProphet and BART networks in the same model average method that NetProphet uses to combine LASSO and DE networks. 

source("~ykang/usr/netprophet_0.1/CODE/modelaverage.r")

args <- commandArgs(trailingOnly=TRUE)
net_np <- as.matrix(read.table(toString(args[1])))
net_bart <- as.matrix(read.table(toString(args[2])))

index <- which(net_np > 0)
net_np[index] <- net_np[index] - min(abs(net_np[index]))
index <- which(net_np < 0)
net_np[index] <- net_np[index] - min(abs(net_np[index]))
net_np <- net_np / max(abs(net_np))

index <- which(net_bart > 0)
net_bart[index] <- net_bart[index] - min(abs(net_bart[index]))
index <- which(net_bart < 0)
net_bart[index] <- net_bart[index] - min(abs(net_bart[index]))
net_bart <- net_bart / max(abs(net_bart))

net_comb <- compute.model.average.new(net_np, net_bart, c(1,1,1,1,1,1,0,0)))

write.table(net_comb, args[length(args)], row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
