args <- commandArgs(trailingOnly=TRUE)
fnPredictedNetwork <- toString(args[1])
fnTfs <- toString(args[2])
fnGenes <- toString(args[3])
fnGoldNetwork <- toString(args[4])
fnOutName <- toString(args[5])
maxRank <- as.integer(args[6])*1000

rank_cutoffs <- seq(maxRank/10, maxRank, maxRank/10)

predictedNetwork <- as.matrix(read.table(fnPredictedNetwork, header=FALSE))
goldNetwork <- as.matrix(read.table(fnGoldNetwork), header=FALSE)
tfs <- as.matrix(read.table(fnTfs))
genes <- as.matrix(read.table(fnGenes))

evaluateEvidSupport <- function(predictedNetwork, goldNetwork, tfs, genes, rank_cutoffs) {
    goldNetwork <- filter_auto(goldNetwork, tfs, genes)
    # predictedNetwork <- filter_auto(predictedNetwork, tfs, genes)

    gold_inds <- which(apply(goldNetwork,1,max)>0)

    res <- matrix(nrow=4, ncol=length(rank_cutoffs), data=0)
    res[1,] <- rep(length(goldNetwork[gold_inds,]), length(rank_cutoffs))
    res[2,] <- rep(sum(goldNetwork[gold_inds,]), length(rank_cutoffs))

    # evaluate all tf edges
    cutoffs <- abs(predictedNetwork)[order(abs(predictedNetwork), decreasing=TRUE)][rank_cutoffs]
	cat(cutoffs, "\n")

    for(j in 1:length(cutoffs)) {
        for (gold_ind in gold_inds) {
            target_inds <- which(abs(predictedNetwork[gold_ind,]) >= cutoffs[j])
            res[3,j] <- res[3,j] + length(target_inds)
            res[4,j] <- res[4,j] + sum(goldNetwork[gold_ind,][target_inds])
        }
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
