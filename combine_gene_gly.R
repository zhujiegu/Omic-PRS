library(magrittr)
library(dplyr)
library(data.table)

# Genetic data with only snps siginificantly related with glycans
load('Intl_snp_gly_training.RData')
SNP_list <- readRDS('SNPs_gly.rds')
dat <- dat[,colnames(dat)%in%SNP_list]
print(dim(dat))
save(dat, gly, file='Intl_snp_gly_sel_training.RData')

load('Intl_snp_test.RData')
dat_t <- dat_t[,colnames(dat_t)%in%SNP_list]
print(dim(dat))
save(dat_t, file='Intl_snp_sel_test.RData')

#####################################################
# # SNP training data
# 
# # Read and manipulate SNP data
# dat <- fread('chr_all_training.raw', header = T, verbose = F)
# 
# psam_training <- read.table('psam_training', header = T)
# rownames(dat) <- psam_training$IID
# 
# snp_id <- colnames(dat)
# 
# snp <- sapply(snp_id, function(e) strsplit(e,'_')[[1]][1])
# colnames(dat) <- snp
# print(dim(dat))
# 
# print(dat[1:6,1:6])
# 
# # Read glycan
# gly <- readRDS("Glycan_sel.rds")
# # gly <- readRDS('/home/z/Data/TwinsUK_CWP/Glycan_sel.rds')
# 
# gly <- gly[match(psam_training$IID, rownames(gly)), ]
# if(!all.equal(rownames(gly), psam_training$IID)) stop('samples do not match')
# 
# print("check point 1")
# gly %<>% scale(scale=F)
# dat %<>% scale(scale = F)
# print("check point 2")
# 
# save(dat, gly, file='Intl_snp_gly_training.RData')
# rm(dat); gc()
# #####################################################
# # SNP test data
# dat_t <- fread('chr_all_test.raw', header = T, verbose = F)
# 
# psam_test <- read.table('psam_test', header = T)
# rownames(dat_t) <- psam_test$IID
# 
# snp_id <- colnames(dat_t)
# 
# snp <- sapply(snp_id, function(e) strsplit(e,'_')[[1]][1])
# colnames(dat_t) <- snp
# print(dim(dat_t))
# 
# print(dat_t[1:6,1:6])
# 
# print("check point 3")
# dat_t %<>% scale(scale = F)
# print("check point 4")
# 
# save(dat_t, file='Intl_snp_test.RData')

#####################################################
# # Summarized gene data
# 
# gene <- readRDS('GeneDat.Rds')
# glycan <- readRDS('Glycan_sel.rds')
# 
# # 2000 training samples
# id <- read.table('id_A_training_cwp.txt', header = T)
# 
# # select samples
# gly <- glycan[match(id$IID, rownames(glycan)),]
# gene <- gene[match(id$IID, rownames(gene)),]
# # check
# dim(gly); dim(gene)
# 
# all.equal(rownames(gly), rownames(gene))
# all.equal(rownames(gly), id$IID)
# 
# gly %<>% scale(scale = F)
# gene %<>% scale(scale = F)
# 
# save(gene, gly, id, file = 'Intl_gene_gly_sel.RData')
# 
# #####################################################
# # gene testing data
# 
# id_test <- read.table('id_AB_test_cwp.txt', header = T)
# gene <- readRDS('GeneDat.Rds')
# gene_test <- gene[match(id_test$IID, rownames(gene)),]
# 
# all.equal(rownames(gene_test), id_test$IID)
# 
# save(gene_test, id_test, file = 'Gene_cwp_test.RData')
# 
# #####################################################
# 
# gene <- readRDS('GeneDat.Rds')
# glycan <- readRDS('Glycan.Rds')
# 
# # overlapping samples
# gly <- glycan[rownames(glycan) %in% rownames(gene),]
# gene <- gene[rownames(gene) %in% rownames(glycan),]
# # check
# dim(gly); dim(gene)
# 
# gly <- gly[order(rownames(gly)),]
# gene <- gene[order(rownames(gene)),]
# 
# all.equal(rownames(gly), rownames(gene))
# 
# gly %<>% scale(scale = F)
# gene %<>% scale(scale = F)
# 
# cwp <- readRDS('cwp.rds')
# all.equal(rownames(gly), cwp$PublicID)
# cwp <- cwp$case
# 
# save(gene, gly, cwp, file = 'Intl_gene_gly.RData')