# Analyze fly network. The analysis scope is the intersection of inferred network and known network compiled from ChIP and PWM datasets respectively. Each column is a cumulative evaluation of the count above certain cutoff of the network edge score.
# 
# Result includes 6 row:
# 1. Number of edges in the intersected network (the size of the network)
# 2. Number of positive edges in the known network (binary data from either ChIP or PWM)
# 3. Number of edges in the inferred network that are above certain cutoff
# 4. Number of supported edges in the inferred network by known network
# 5. Number of supported TFs that have >= 1 supported edges


args <- commandArgs(trailingOnly=TRUE)
fnPredictedNetwork <- toString(args[1])
fnTfs <- toString(args[2])
fnGenes <- toString(args[3])
fnGoldNetwork <- toString(args[4])
fnOutName <- toString(args[5])
maxRank <- as.integer(args[6])*1000
if (length(args) < 7) {
    numBins <- 10
} else {
    numBins <- as.integer(args[7])
}

cat("Max Ranking:", maxRank, ", Number of bins:", numBins,"\n")

rank_cutoffs <- seq(maxRank/numBins, maxRank, maxRank/numBins)

predictedNetwork <- as.matrix(read.table(fnPredictedNetwork, header=FALSE))
goldNetwork <- as.matrix(read.table(fnGoldNetwork), header=FALSE)
tfs <- as.matrix(read.table(fnTfs))
genes <- as.matrix(read.table(fnGenes))

evaluateEvidSupport <- function(predictedNetwork, goldNetwork, tfs, genes, rank_cutoffs) {
    goldNetwork <- filter_auto(goldNetwork, tfs, genes)
    # predictedNetwork <- filter_auto(predictedNetwork, tfs, genes)

    gold_inds <- which(apply(goldNetwork,1,max)>0)

    res <- matrix(nrow=5, ncol=length(rank_cutoffs), data=0)
    res[1,] <- rep(length(goldNetwork[gold_inds,]), length(rank_cutoffs))
    res[2,] <- rep(sum(goldNetwork[gold_inds,]), length(rank_cutoffs))

    # evaluate all tf edges
    cutoffs <- abs(predictedNetwork)[order(abs(predictedNetwork), decreasing=TRUE)][rank_cutoffs]
	cat("Edge cutoffs:", cutoffs, "\n")

    for(j in 1:length(cutoffs)) {
        tfs_supported <- c()
        for (gold_ind in gold_inds) {
            target_inds <- which(abs(predictedNetwork[gold_ind,]) >= cutoffs[j])
            if (length(target_inds) > 0) {
                tfs_supported <- c(tfs_supported, tfs[gold_ind])
            }
            res[3,j] <- res[3,j] + length(target_inds)
            res[4,j] <- res[4,j] + sum(goldNetwork[gold_ind,][target_inds])
        }
        res[5,j] <- length(tfs_supported)
        
        # predictedSubNet <- abs(predictedNetwork)[gold_inds,]
        # inds <- which(predictedSubNet > cutoffs[j])
        # res[3,j] <- length(inds)
        # res[4,j] <- sum(goldNetwork[inds])
    }

    # # evaluate chip tf edges
    # cutoffs <- abs(predictedNetwork[gold_inds,])[order(abs(predictedNetwork[gold_inds,]), decreasing=TRUE)][rank_cutoffs]
    # for(j in 1:length(cutoffs)) {
    #     if(j == 1) {
    #             inds <- which(abs(predictedNetwork)[gold_inds,] >= cutoffs[j])
    #     } else {
    #             inds <- which(abs(predictedNetwork)[gold_inds,] >= cutoffs[j] & abs(predictedNetwork)[gold_inds,] < cutoffs[j-1])
    #     }
    #     res[3,j] <- sum(goldNetwork[gold_inds,][inds])/length(inds)
    # }

    write.table(res, file=fnOutName, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
}

filter_auto <- function(network, tfs, genes){
    select <- cbind(1:length(tfs), match(tfs,genes))
    select <- select[which(!is.na(select[,2])),]
    network[select] <- 0
    network
}

evaluateEvidSupport(predictedNetwork, goldNetwork, tfs, genes, rank_cutoffs)
