######################################################################
############ GWAS ########

#### get EUR ancestry specific PCs

sed 's/_/ /g' plink2.king.cutoff.in.id > plink2.king.cutoff.in.id_edit

/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/software/plink2 --bfile /sc/arion/projects/psychgen/cotea02_prset/chronicpain_mr/biome_chronicpain_gwas/data/ancestry/1000G_BioMe_forFlashPCA --keep plink2.king.cutoff.in.id_edit --make-bed --out /sc/arion/projects/psychgen/cotea02_prset/chronicpain_mr/biome_chronicpain_gwas/data/ancestry/BioMe_EUR_forFlashPCA

ml flashpca
flashpca --bfile /sc/arion/projects/psychgen/cotea02_prset/chronicpain_mr/biome_chronicpain_gwas/data/ancestry/BioMe_EUR_forFlashPCA --ndim 50


#### build phenofile (PCs, IID, Dx, Age, Sex, Chip)

sed 's/0x//g' chronic_pain_broad_patientids.txt > tmp && mv tmp chronic_pain_broad_patientids.txt
sed 's/0x//g' chronic_pain_narrow_patientids.txt > tmp && mv tmp chronic_pain_narrow_patientids.txt

library(data.table)
library(tidyverse)
setwd('/sc/arion/projects/psychgen/cotea02_prset/chronicpain_mr/biome_chronicpain_gwas/data/')
broad <- read.table('./biome_phenotype/chronic_pain_broad_patientids.txt',fill=T)
narrow <- read.table('./biome_phenotype/chronic_pain_narrow_patientids.txt',fill=T)
ids <- read.table('./ancestry/biome_eur_ancestry_ids.txt')
pcs <- read.table('./ancestry/pcs.txt',header=T)
dems <- fread('/sc/arion/projects/data-ark/CBIPM/datamart/BioMe_Anonymized/2023-04-03/person.txt',sep="|") %>% as.data.frame()

pcs.eur <- pcs[pcs$IID %in% ids$V1,]
pcs.eur$broad <- ifelse(pcs.eur$IID %in% broad$V1,"2","1")
pcs.eur$narrow <- ifelse(pcs.eur$IID %in% narrow$V1,"2","1")

dems$IID <- gsub("0x","",dems$person_id)
dems$Age <- 2025-dems$year_of_birth
dems <- dems[,colnames(dems) %in% c("IID","Age")]
temp <- merge(pcs.eur,dems,by="IID")

fam <- read.table('/sc/arion/projects/data-ark/CBIPM/Microarray/combined/genotyped_V2/GSA_GDA_HG37_genotyped.fam')
fam$IID <- fam$V1
fam$Sex <- fam$V5
fam$CHIP <- fam$V7
fam <- fam[,colnames(fam) %in% c('Sex','IID',"CHIP")]

master <- merge(temp,fam,by="IID")
write.table(master,"./biome_phenotype/biome_phenofile_forgwas.txt",quote=F,row.names=F)


### reformat phenofile IIDs

awk '{print $1"_"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' ./data/biome_phenotype/biome_phenofile_forgwas.txt > tmp && mv tmp ./data/biome_phenotype/biome_phenofile_forgwas.txt


#### run GWAS

path="/sc/arion/projects/psychgen/cotea02_prset/chronicpain_mr/biome_chronicpain_gwas"

/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/software/plink2 \
--pfile ${path}/data/geno/GSA_GDA_TOPMED_allchr.8bit_EUR_foundersonly \
--glm hide-covar \
--covar ${path}/data/biome_phenotype/biome_phenofile_forgwas.txt \
--covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Age,Sex,CHIP \
--pheno ${path}/data/biome_phenotype/biome_phenofile_forgwas.txt \
--pheno-name broad \
--covar-variance-standardize \
--out ${path}/results/biome_eur_broadchronicpain_gwas



path="/sc/arion/projects/psychgen/cotea02_prset/chronicpain_mr/biome_chronicpain_gwas"

/sc/arion/projects/psychgen/cotea02_prset/geneoverlap/software/plink2 \
--pfile ${path}/data/geno/GSA_GDA_TOPMED_allchr.8bit_EUR_foundersonly \
--glm hide-covar \
--covar ${path}/data/biome_phenotype/biome_phenofile_forgwas.txt \
--covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Age,Sex,CHIP \
--pheno ${path}/data/biome_phenotype/biome_phenofile_forgwas.txt \
--pheno-name narrow \
--covar-variance-standardize \
--out ${path}/results/biome_eur_narrowchronicpain_gwas
