args <- commandArgs(trailingOnly=TRUE)
fnPredictedNetwork <- toString(args[1])
fnTfs <- toString(args[2])
fnGenes <- toString(args[3])
fnGoldNetwork <- toString(args[4])
fnOutName <- toString(args[5])

# rank_cutoffs <- seq(2000, 20000, 2000)
rank_cutoffs <- seq(4000, 40000, 4000)
# rank_cutoffs <- seq(20000, 200000, 20000)

predictedNetwork <- as.matrix(read.table(fnPredictedNetwork, header=TRUE))
goldNetwork <- as.matrix(read.table(fnGoldNetwork))
tfs <- as.matrix(read.table(fnTfs))
genes <- as.matrix(read.table(fnGenes))

evaluateEvidSupport <- function(predictedNetwork, goldNetwork, tfs, genes, rank_cutoffs) {
        goldNetwork <- filter_auto(goldNetwork, tfs, genes)
        # predictedNetwork <- filter_auto(predictedNetwork, tfs, genes)

        gold_inds <- which(apply(goldNetwork,1,max)>0)
        # cat(gold_inds, "\n")
        rand <- sum(goldNetwork[gold_inds,]) / length(goldNetwork[gold_inds,])
        res <- matrix(nrow=3, ncol=length(rank_cutoffs), data=0)
        res[1,] <- rep(rand, length(rank_cutoffs))

        # evaluate all tf edges
        cutoffs <- abs(predictedNetwork)[order(abs(predictedNetwork), decreasing=TRUE)][rank_cutoffs]
        for(j in 1:length(cutoffs)) {
                if(j == 1) {
                        inds <- which(abs(predictedNetwork)[gold_inds,] >= cutoffs[j])
                } else {
                        inds <- which(abs(predictedNetwork)[gold_inds,] >= cutoffs[j] & abs(predictedNetwork)[gold_inds,] < cutoffs[j-1])
                }
                res[2,j] <- sum(goldNetwork[gold_inds,][inds])/length(inds)
        }

        # evaluate chip tf edges
        cutoffs <- abs(predictedNetwork[gold_inds,])[order(abs(predictedNetwork[gold_inds,]), decreasing=TRUE)][rank_cutoffs]
        for(j in 1:length(cutoffs)) {
                if(j == 1) {
                        inds <- which(abs(predictedNetwork)[gold_inds,] >= cutoffs[j])
                } else {
                        inds <- which(abs(predictedNetwork)[gold_inds,] >= cutoffs[j] & abs(predictedNetwork)[gold_inds,] < cutoffs[j-1])
                }
                res[3,j] <- sum(goldNetwork[gold_inds,][inds])/length(inds)
        }

        write.table(res, file=fnOutName, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
}

filter_auto <- function(network, tfs, genes){
        select <- cbind(1:length(tfs), match(tfs,genes))
        select <- select[which(!is.na(select[,2])),]
        network[select] <- 0
        network
}

evaluateEvidSupport(predictedNetwork, goldNetwork, tfs, genes, rank_cutoffs)

# evaluateEvidSupport(predictedNetwork, goldNetwork[match(intersect(combine_tfs, tfs), tfs), match(intersect(combine_genes, genes), genes)], combine_tfs, combine_genes, seq(2000, 20000, 2000))
