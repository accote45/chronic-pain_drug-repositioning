
library(tidyverse)
library(data.table)
library(plyr)
library(dplyr)
library(EMUDRA)
library(data.table)

source("/sc/arion/projects/mscic1/results/sigmap_prelim/emudra/Source/EMUDRA/R/eXtremeLogFC.R")
source("/sc/arion/projects/mscic1/results/sigmap_prelim/emudra/Source/EMUDRA/R/eXScores.R")
source("/sc/arion/projects/mscic1/results/sigmap_prelim/emudra/Source/EMUDRA/R/XCor.R")
source("/sc/arion/projects/mscic1/results/sigmap_prelim/emudra/Source/EMUDRA/R/XCorrelation.R")
source("/sc/arion/projects/mscic1/results/sigmap_prelim/emudra/Source/EMUDRA/R/XCosine.R")
source("/sc/arion/projects/mscic1/results/sigmap_prelim/emudra/Source/EMUDRA/R/XCos.R")
source("/sc/arion/projects/mscic1/results/sigmap_prelim/emudra/Source/EMUDRA/R/XSum.R")

cmap_score <- function(sig_up, sig_down, drug_signature) {

        #the old function does not support the input list with either all up genes or all down genes, this new function attempts to addess this, not fully validated
        num_genes <- nrow(drug_signature)
        ks_up <- 0
        ks_down <- 0
        connectivity_score <- 0

        drug_signature[,"rank"] <- rank(-drug_signature[,2])

        up_tags_rank <- merge(drug_signature, sig_up, by.x = "ids", by.y = 1)
        down_tags_rank <- merge(drug_signature, sig_down, by.x = "ids", by.y = 1)

        up_tags_position <- sort(up_tags_rank$rank)
        down_tags_position <- sort(down_tags_rank$rank)

        num_tags_up <- length(up_tags_position)
        num_tags_down <- length(down_tags_position)

        # 
        if(num_tags_up > 1) {
                a_up <- 0
                b_up <- 0

                a_up <- max(sapply(1:num_tags_up,function(j) {
                        j/num_tags_up - up_tags_position[j]/num_genes
                }))
                b_up <- max(sapply(1:num_tags_up,function(j) {
                        up_tags_position[j]/num_genes - (j-1)/num_tags_up
                }))

                if(a_up > b_up) {
                        ks_up <- a_up
                } else {
                        ks_up <- -b_up
                }
        }else{
                ks_up <- 0
        }

        if (num_tags_down > 1){

                a_down <- 0
                b_down <- 0

                a_down <- max(sapply(1:num_tags_down,function(j) {
                        j/num_tags_down - down_tags_position[j]/num_genes
                }))
                b_down <- max(sapply(1:num_tags_down,function(j) {
                        down_tags_position[j]/num_genes - (j-1)/num_tags_down
                }))

                if(a_down > b_down) {
                        ks_down <- a_down
                } else {
                        ks_down <- -b_down
                }
        }else{
                ks_down <- 0
        }

        if (ks_up == 0 & ks_down != 0){ #only down gene inputed
                connectivity_score <- -ks_down
        }else if (ks_up !=0 & ks_down == 0){ #only up gene inputed
                connectivity_score <- ks_up
        }else if (sum(sign(c(ks_down,ks_up))) == 0) {
                connectivity_score <- ks_up - ks_down # different signs
        }

        return(connectivity_score)
}


args <- commandArgs(trailingOnly=TRUE)
number=args[1]
perm=args[2]


mat <- read.table(paste0('/sc/arion/projects/mscic1/results/sigmap_prelim/cmap_lincs_2020/level5_beta_trt_cp_n720216x12328_mat.',number,'.gz'),header=T)

ids <- read.table('/sc/arion/projects/mscic1/results/sigmap_prelim/sigmap_manual_results/chronic_pain/permuted_scores/mat.ids.txt',sep=" ",header=T)

ids$gene <- as.character(as.numeric(ids$gene))

rows <- read.table("/sc/arion/projects/mscic1/results/sigmap_prelim/cmap_lincs_2020/rdesc",header=TRUE)

mat$'gene' <- rows$id
mat$gene <- as.character(as.numeric(mat$gene))
mat <- mat %>% left_join(ids, by=c("gene" = "gene"))
mat <- mat[!duplicated(mat$hgnc_symbol),]
mat <- na.omit(mat)
mat <- as.data.frame(mat)


############################################################################
### shuffle disease z-scores, compare to drug signatures

##TWAS QVAL
rand.qval <- function(fdr){
	twas <- read.table('/sc/arion/projects/mscic1/results/sigmap_prelim/s_predixcan/pain_johnston/gtexv8_Muscle_Skeletal.csv_with_qval.txt',header=T)
	twas <- twas[!is.na(twas$zscore),]
	twas <- transform(twas,gene_rand=sample(gene))

	twas.qval <- read.table(paste0('/sc/arion/projects/mscic1/results/sigmap_prelim/s_predixcan/pain_johnston/gtexv8_Muscle_Skeletal_',fdr,'qval.txt'),header=F)
	colnames(twas.qval) <- c('gene','gene_name','zscore','effect_size','pvalue','VAR_g','pred_perf_R2','pred_perf_p','pred_perf_q','n_snps_used','n_snps_in_cov','n_snps_in_model','qval')
	twas.qval$up_down[as.numeric(as.character(twas.qval$zscore))>0] <- "up"
	twas.qval$up_down[as.numeric(as.character(twas.qval$zscore))<0] <- "down"
	dz_signature <- twas.qval
	dz_genes_up <- subset(dz_signature,up_down=="up",select=c("gene_name","zscore"))
	dz_genes_down <- subset(dz_signature,up_down=="down",select=c("gene_name","zscore"))

	gene_list_twas <- subset(twas,select=c("gene_rand","zscore"))
	gene_list_twas_up <- gene_list_twas[gene_list_twas$zscore>0,]
	gene_list_twas_down <- gene_list_twas[gene_list_twas$zscore<0,]


	#random_input_signature_genes <- sample(gene_list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
	rand_dz_gene_up <- gene_list_twas_up[1:nrow(dz_genes_up),]
	rand_dz_gene_down <- gene_list_twas_down[1:nrow(dz_genes_down),]
	mat.x <- mat
	rownames(mat.x) <- mat.x$ensembl_gene_id
	mat.x$gene <- NULL
	mat.x$hgnc_symbol <- NULL
	mat.x$ensembl_gene_id <- NULL

######## eXtreme scores
	x.scores.qval.dat <- eXScores(mat.x, rand_dz_gene_up, rand_dz_gene_down)
	colnames(x.scores.qval.dat) <- c("xsum","xcos","xcor","xspe")
	x.scores.qval.dat <- as.data.frame(x.scores.qval.dat) %>% rownames_to_column('signature')

	qval.long <- gather(x.scores.qval.dat,score_type,value,xsum:xspe)
	qval.long$twas_sig <- fdr
	assign(paste("qval.long",fdr,sep=""),qval.long,envir=parent.frame())

######## WCS

	gene_list <- rownames(mat.x)

	disease_cmap_scores <- sapply(1:ncol(mat.x),function(exp_id) {
	print(paste("Computing score for disease against cmap_experiment_id =",exp_id))
	cmap_exp_signature <- cbind(gene_list,subset(mat.x,select=exp_id))
	colnames(cmap_exp_signature) <- c("ids","rank")
	cmap_score(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature)
	},simplify=F)

	wcs <- as.data.frame(unlist(disease_cmap_scores))
	wcs$rank <- rank(wcs[,1])
	colnames(wcs) <- c("score","rank")
	wcs$perturbation <- colnames(mat.x)
	wcs <- wcs[order(wcs$rank),]

	wcs$twas_sig <- num
	assign(paste("wcs_",num,sep=""),wcs,envir=parent.frame())
}

seq <- c('0.05','0.1','0.2')
for (num in seq){
rand.qval(num)
}


### TWAS finemapped
twas <- read.table('/sc/arion/projects/mscic1/results/sigmap_prelim/s_predixcan/finemap_focus/chronic_pain_singletissue_reference/pain_johnston/ChronicPain_Muscle_Skeletal_credsets.tsv',sep="\t",header=T)
twas <- twas[!(twas$block=="block"),]
twas$twas_z_pop1 <- as.numeric(as.character(twas$twas_z_pop1))
twas$up_down[twas$twas_z_pop1>0] <- "up"
twas$up_down[twas$twas_z_pop1<0] <- "down"
dz_signature <- twas

twas <- transform(twas,zscore_rand=sample(twas_z_pop1))

dz_genes_up <- subset(dz_signature,up_down=="up",select=c("mol_name"))
dz_genes_down <- subset(dz_signature,up_down=="down",select=c("mol_name"))

	gene_list_twas <- subset(twas,select=c("ens_gene_id","zscore_rand"))

	#random_input_signature_genes <- sample(gene_list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
	rand_dz_gene_up <- gene_list_twas[1:nrow(dz_genes_up),]
	rand_dz_gene_down <- gene_list_twas[(nrow(dz_genes_up)+1):nrow(gene_list_twas),]
	mat.x <- mat
	rownames(mat.x) <- mat.x$ensembl_gene_id
	mat.x$gene <- NULL
	mat.x$hgnc_symbol <- NULL
	mat.x$ensembl_gene_id <- NULL

######## eXtreme scores
	x.scores.qval.dat <- eXScores(mat.x, rand_dz_gene_up, rand_dz_gene_down)
	colnames(x.scores.qval.dat) <- c("xsum","xcos","xcor","xspe")
	x.scores.qval.dat <- as.data.frame(x.scores.qval.dat) %>% rownames_to_column('signature')

	qval.long <- gather(x.scores.qval.dat,score_type,value,xsum:xspe)
	qval.long$twas_sig <- "credset"
	assign(paste("credset.long"),qval.long,envir=parent.frame())

######## WCS

	gene_list <- rownames(mat.x)

	disease_cmap_scores <- sapply(1:ncol(mat.x),function(exp_id) {
	print(paste("Computing score for disease against cmap_experiment_id =",exp_id))
	cmap_exp_signature <- cbind(gene_list,subset(mat.x,select=exp_id))
	colnames(cmap_exp_signature) <- c("ids","rank")
	cmap_score(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature)
	},simplify=F)

	wcs <- as.data.frame(unlist(disease_cmap_scores))
	wcs$rank <- rank(wcs[,1])
	colnames(wcs) <- c("score","rank")
	wcs$perturbation <- colnames(mat.x)
	wcs <- wcs[order(wcs$rank),]

	wcs$twas_sig <- "credset"
	assign(paste("wcs_credset"),wcs,envir=parent.frame())


master <- bind_rows(wcs_0.05,wcs_0.1,wcs_0.2,wcs_credset)
write.table(master,paste0("/sc/arion/projects/mscic1/results/sigmap_prelim/sigmap_manual_results/chronic_pain/permuted_scores_onescramble/pain_khoury/skeletal_muscle/wcs_perm_",number,"_",perm,".txt"),quote=FALSE,row.names=FALSE)


master2 <- bind_rows(qval.long0.05,qval.long0.1,qval.long0.2,credset.long)
write.table(master2,paste0("/sc/arion/projects/mscic1/results/sigmap_prelim/sigmap_manual_results/chronic_pain/permuted_scores_onescramble/pain_khoury/skeletal_muscle/extreme_perm_",number,"_",perm,".txt"),quote=F,row.names=F)




