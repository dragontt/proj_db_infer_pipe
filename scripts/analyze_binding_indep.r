args <- commandArgs(trailingOnly = TRUE)
network <- toString(args[1])
chip_net <- toString(args[2])
pwm_net <- toString(args[3])

net <- read.table(network,header=TRUE)
pdna.inter <- read.table(chip_net,header=TRUE)
bsinfo <- read.table(pwm_net,header=TRUE)

source("~bhaynes/projects/network_inference/util/util.r")

chip.bp.np.setsizes <- c()

orf_universe <- intersect(union(unique(bsinfo$TARGET),unique(pdna.inter$TARGET)),unique(net$TARGET))
reg_pdna.inter <- intersect(unique(as.character(pdna.inter$REGULATOR)),unique(as.character(net$REGULATOR)))
reg_bsinfo <- intersect(intersect(unique(as.character(bsinfo$REGULATOR)),unique(as.character(net$REGULATOR))), unique(as.character(pdna.inter$REGULATOR)))
interaction_pdna.inter_cnt <- length(reg_pdna.inter) * length(orf_universe)
interaction_bsinfo_cnt <- length(reg_bsinfo) * length(orf_universe)

# np_conf_cutoffs <- net$CONFIDENCE[order(net$CONFIDENCE,decreasing=TRUE)][seq(2000,20000,by=2000)]
np_conf_cutoffs <- net$CONFIDENCE[order(net$CONFIDENCE,decreasing=TRUE)][seq(4000,40000,by=4000)]
# np_conf_cutoffs <- net$CONFIDENCE[order(net$CONFIDENCE,decreasing=TRUE)][seq(20000,200000,by=20000)]
for (np_conf_cutoff in np_conf_cutoffs)
{
	for (chip_desired_recall in c(0.1)) # c(0.1,0.33,0.5)
	{
		bp.cutoff <- c()
		chip.bp.np.sets <- c()

		i <- 1
		for (regulator in reg_pdna.inter)
		{
			# targets <- union(as.character(pdna.inter$TARGET[which(pdna.inter$REGULATOR==regulator)]),as.character(bsinfo$TARGET[which(bsinfo$REGULATOR==regulator)]))
			targets <- orf_universe
			pdna.evid <- rep(0,times=length(targets))
			binding.evid.max <- rep(0,times=length(targets))
			binding.evid.sum <- rep(0,times=length(targets))
			pdna.evid[match(as.character(pdna.inter$TARGET[which(pdna.inter$REGULATOR==regulator)]),targets)] <- 1
			binding.evid.max[match(as.character(bsinfo$TARGET[which(bsinfo$REGULATOR==regulator)]),targets)] <- bsinfo$MAXP[which(bsinfo$REGULATOR==regulator)]
			binding.evid.sum[match(as.character(bsinfo$TARGET[which(bsinfo$REGULATOR==regulator)]),targets)] <- bsinfo$SUMP[which(bsinfo$REGULATOR==regulator)]
			cat(regulator,length(targets),"\n")
			binding.evid.comb <- apply(rbind(binding.evid.max,binding.evid.sum),2,max)	
			prc.comb <- compute.prc(binding.evid.comb,pdna.evid)
		
			if 	(length(which(prc.comb[,1]<=chip_desired_recall))>0)
			{
				bp.cutoff[i] <- binding.evid.comb[sort.list(binding.evid.comb,decreasing=TRUE)][length(which(prc.comb[,1]<=chip_desired_recall))]
			}
			else
			{
				bp.cutoff[i] <- 1
			}
		
			predicted.targets <- as.character(net$TARGET[which(net$REGULATOR==regulator)])[which(net$CONFIDENCE[which(net$REGULATOR==regulator)]>np_conf_cutoff)]
			predicted.targets <- intersect(predicted.targets,orf_universe) # Narrow to space under consideration

			bp.targets <- targets[which(binding.evid.comb>bp.cutoff[i])]
			pdna.targets <- targets[which(pdna.evid>0)]
			
			chip.bp.np.sets <- rbind(chip.bp.np.sets,c(length(bp.targets),sum(pdna.evid),length(predicted.targets),length(intersect(bp.targets,pdna.targets)),length(intersect(bp.targets,predicted.targets)),length(intersect(predicted.targets,pdna.targets)),length(intersect(bp.targets,intersect(predicted.targets,pdna.targets)))))
	
			i <- i + 1
		}
		
		chip.bp.np.setsize <- apply(chip.bp.np.sets,2,sum)
		chip.bp.np.setsizes <- rbind(chip.bp.np.setsizes,chip.bp.np.setsize)
	}
}

# filename <- paste("chip.bp.np.set.sizes.top2to20k.",strsplit(network,"[./]")[[1]][4],".txt", sep="")
filename <- paste("chip.bp.np.set.sizes.top4to40k.",strsplit(network,"[./]")[[1]][4],".txt", sep="")
# filename <- paste("chip.bp.np.set.sizes.top20to200k.",strsplit(network,"[./]")[[1]][4],".txt", sep="")
write.table(cbind(chip.bp.np.setsizes, rep(interaction_bsinfo_cnt, length=nrow(chip.bp.np.setsizes)), rep(interaction_pdna.inter_cnt, length=nrow(chip.bp.np.setsizes))),file=filename,col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')

