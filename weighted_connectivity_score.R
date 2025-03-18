

## partly adapted from https://github.com/Bin-Chen-Lab/HCC_NEN/

## script for example tissue (GTEx-ACC)

############# Sig Map Calculation ############

library(data.table)
library(tidyverse)
library(plyr)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
number=args[1]

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


# Load in the signature for the subset_comparison_id
mat <- fread(paste0('/sc/arion/projects/psychgen/alanna/sigmap_prelim/cmap_lincs_2020/level5_beta_trt_cp_n720216x12328_mat.',number,'.gz'))

mat <- as.data.frame(mat)
rows <- read.table("/sc/arion/projects/psychgen/alanna/sigmap_prelim/cmap_lincs_2020/rdesc",header=TRUE)

library(biomaRt)
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",mirror="uswest")
Ensemble2HGNC <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                       filters = "entrezgene_id", values = rows$id,
                       mart = mart, uniqueRows = TRUE)
Ensemble2HGNC <- Ensemble2HGNC[complete.cases(Ensemble2HGNC),]
Ensemble2HGNC <- Ensemble2HGNC[!duplicated(Ensemble2HGNC$hgnc_symbol),]

Ensemble2HGNC$gene <- Ensemble2HGNC$hgnc_symbol
Ensemble2HGNC$entrezgene_id <- as.character(as.numeric(Ensemble2HGNC$entrezgene_id))

rownames(mat) <- rows$id
mat <- mat %>% rownames_to_column('gene')
mat <- mat %>% left_join(Ensemble2HGNC, by=c("gene" = "entrezgene_id"))
mat <- mat[,-c(1,(ncol(mat)-1))]

mat <- mat%>% drop_na(gene.y)
# 11,349 genes included in cmap resource after converting from entrez to hgnc id

################################################
######## TWAS QVAL #############################
seq <- c('0.01','0.05','0.1','0.2')

for (num in seq){
    twas <- read.table(paste0('/sc/arion/projects/psychgen/alanna/sigmap_prelim/s_predixcan/chronic_pain/gtexv8_Brain_Anterior_cingulate_cortex_BA24_',num,'qval.txt'),header=F)
    colnames(twas) <- c('gene','gene_name','zscore','effect_size','pvalue','VAR_g','pred_perf_R2','pred_perf_p','pred_perf_q','n_snps_used','n_snps_in_cov','n_snps_in_model','qval')
    twas$up_down[as.numeric(as.character(twas$zscore))>0] <- "up"
    twas$up_down[as.numeric(as.character(twas$zscore))<0] <- "down"
    dz_signature <- twas
    dz_genes_up <- subset(dz_signature,up_down=="up",select=c("gene_name"))
    dz_genes_down <- subset(dz_signature,up_down=="down",select=c("gene_name"))

    gene_list <- subset(mat,select=ncol(mat))
    cmap_signatures <- mat[,1:(ncol(mat)-1)] 

    colnames(dz_genes_up) <- "GeneID"
    colnames(dz_genes_down) <- "GeneID"

    disease_cmap_scores <- sapply(1:ncol(cmap_signatures),function(exp_id) {
    print(paste("Computing score for disease against cmap_experiment_id =",exp_id))
    cmap_exp_signature <- cbind(gene_list,subset(cmap_signatures,select=exp_id))
    colnames(cmap_exp_signature) <- c("ids","rank")
    cmap_score(dz_genes_up,dz_genes_down,cmap_exp_signature)
    },simplify=F)

    wcs <- as.data.frame(unlist(disease_cmap_scores))
    wcs$rank <- rank(wcs[,1])
    colnames(wcs) <- c("score","rank")
    wcs$perturbation <- colnames(cmap_signatures)
    wcs <- wcs[order(wcs$rank),]

    wcs$twas_sig <- num
    assign(paste("wcs_",num,sep=""),wcs)
}


################################################
######## TWAS FINEMAPPED ########################
twas <- read.table('/sc/arion/projects/psychgen/alanna/sigmap_prelim/s_predixcan/finemap_focus/chronic_pain_singletissue_reference/ChronicPain_Brain_Anterior_cingulate_cortex_BA24_credsets.tsv',sep="\t",header=T)
twas <- twas[!(twas$block=="block"),]
twas$twas_z_pop1 <- as.numeric(as.character(twas$twas_z_pop1))
twas$up_down[twas$twas_z_pop1>0] <- "up"
twas$up_down[twas$twas_z_pop1<0] <- "down"
dz_signature <- twas
dz_genes_up <- subset(dz_signature,up_down=="up",select=c("mol_name"))
dz_genes_down <- subset(dz_signature,up_down=="down",select=c("mol_name"))

gene_list <- subset(mat,select=ncol(mat))
    cmap_signatures <- mat[,1:(ncol(mat)-1)] 

    colnames(dz_genes_up) <- "GeneID"
    colnames(dz_genes_down) <- "GeneID"

    disease_cmap_scores <- sapply(1:ncol(cmap_signatures),function(exp_id) {
    print(paste("Computing score for disease against cmap_experiment_id =",exp_id))
    cmap_exp_signature <- cbind(gene_list,subset(cmap_signatures,select=exp_id))
    colnames(cmap_exp_signature) <- c("ids","rank")
    cmap_score(dz_genes_up,dz_genes_down,cmap_exp_signature)
    },simplify=F)

    wcs <- as.data.frame(unlist(disease_cmap_scores))
    wcs$rank <- rank(wcs[,1])
    colnames(wcs) <- c("score","rank")
    wcs$perturbation <- colnames(cmap_signatures)
    wcs <- wcs[order(wcs$rank),]

    wcs$twas_sig <- "credset"
    wcs_credset <- wcs

master <- bind_rows(wcs_0.01,wcs_0.05,wcs_0.1,wcs_0.2,wcs_credset)

#save(disease_cmap_scores,file=paste('chronicpain_psychencode_wcs.RData',sep=""))
write.table(master,paste0("/sc/arion/projects/psychgen/alanna/sigmap_prelim/sigmap_manual_results/chronic_pain/wcs/chronicpain_acc_wcs_dat_",number,".txt"),quote=FALSE,row.names=FALSE)


