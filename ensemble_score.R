
library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)

############################################################################
### read in files, subset for hiq signatures
ids <- read.table('/sc/arion/projects/mscic1/results/sigmap_prelim/cmap_lincs_2020/hiq_signature_ids.txt',header=T)

extreme_scores <- list()
for (num in c(seq(1,4))){
extreme_scores[[num]] <- read.table(paste0('/sc/arion/projects/mscic1/results/sigmap_prelim/sigmap_manual_results/chronic_pain/extreme_scores/pain_khoury/dgn_',num,".txt"),header=TRUE)
extreme_scores[[num]] <- extreme_scores[[num]][extreme_scores[[num]]$twas_sig!="0.01",]
extreme_scores[[num]] <- extreme_scores[[num]][extreme_scores[[num]]$signature %in% ids$sig_id,]
}

extreme_scores_master <- ldply(extreme_scores,rbind)

wcs <- read.table('/sc/arion/projects/mscic1/results/sigmap_prelim/sigmap_manual_results/chronic_pain/wcs/pain_khoury/chronicpain_dgn_wcs_wcs_master_fin',header=TRUE,sep="\t")
wcs <- wcs[wcs$twas_sig!="0.01",]
### removed FDR 1% connectivity scores (too few genes tested)

######################
### extreme scores
### normalize each score
extreme_scores_master$id <- paste0(extreme_scores_master$score_type,"_",extreme_scores_master$twas_sig)

norm <- function(x) {
    (x - mean.val) / sd.val
  }

norm.master <- data.frame()
for (i in unique(extreme_scores_master$id)){
	sub <- extreme_scores_master[extreme_scores_master$id==i,]
	mean.val <- mean(sub$value)
	sd.val <- sd(sub$value)
	sub$norm_score <- sapply(sub$value, norm)
	norm.master <- rbind(norm.master,sub)
	#assign(paste0("norm_dat_",i),sub)
}

### average normalized scores across TWAS significance thresholds

avg.norm.master.extreme <- data.frame()
for (i in unique(norm.master$score_type)){
	sub <- norm.master[norm.master$score_type==i,]
	avg.norm <- aggregate(norm_score~signature,sub,FUN=mean)
	avg.norm$score_type <- i
	avg.norm.master.extreme <- rbind(avg.norm.master.extreme,avg.norm)
}


######################

norm <- function(x) {
    (x - mean.val) / sd.val
  }

norm.master.wcs <- data.frame()
for (i in unique(wcs$twas_sig)){
	sub <- wcs[wcs$twas_sig==i,]
	mean.val <- mean(sub$score)
	sd.val <- sd(sub$score)
	sub$norm_score <- sapply(sub$score, norm)
	norm.master.wcs <- rbind(norm.master.wcs,sub)
	#assign(paste0("norm_dat_",i),sub)
}

### average normalized scores across TWAS significance thresholds

avg.norm.master.wcs <- aggregate(norm_score~perturbation,norm.master.wcs,FUN=mean)
avg.norm.master.wcs$score_type <- "wcs"
colnames(avg.norm.master.wcs) <- c('signature','norm_score','score_type')

avg.norm.master <- rbind(avg.norm.master.wcs,avg.norm.master.extreme)



### calculate ensemble score, average across scores
ensemble.master <- aggregate(norm_score~signature,avg.norm.master,FUN=mean)

map <- read.table('/sc/arion/projects/mscic1/results/sigmap_prelim/cmap_lincs_2020/siginfo_beta_hiq_only.txt',sep="\t",header=TRUE)
map <- subset(map,select=c("cmap_name","sig_id","pert_id","cell_iname"))
map$signature <- map$sig_id

fin <- merge(ensemble.master,map,by="signature")
fin <- fin[order(fin$norm_score),]
fin$rank <- rank(fin$norm_score,ties.method="random")
fin$sig_id <- NULL

write.table(fin,"/sc/arion/projects/mscic1/results/sigmap_prelim/sigmap_manual_results/chronic_pain/ensemble_scores/pain_khoury/dgn_ensemble.txt",quote=F,row.names=F,sep="\t")





