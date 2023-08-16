library(OmicsPLS)
library(dplyr)

id_train <- read.table('id_train')
id_test <- read.table('id_test')
colnames(id_train) <- colnames(id_test) <- 'IID'

# function to calculate joint scores in new samples for O2PLS
o2_Tt <- function(fit, dat){
  To <- dat %*% fit$W_Yosc
  dat <- dat - To %*% t(fit$P_Yosc.)
  return(dat %*% fit$W.)
}
###########################################################################################
# # unscreened
# mtb <- readRDS('mtb_sel_bmi.rds')
# snp <- readRDS("../genotypes/pgen/dat_snp_orcades.rds")
# gene <- readRDS("../genotypes/pgen/GeneData_orcades.rds")
# 
# snp <- scale(snp, scale = T)
# gene <- scale(gene, scale = T)
# 
# # test data
# gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
# snp_t <- snp[rownames(snp)%in% id_test$IID, ] %>% scale(scale=F)
# dim(gene_t); dim(snp_t)
# 
# # training data
# mtb <- mtb %>% filter(IID %in% id_train$IID)
# gene <- gene[rownames(gene)%in% mtb$IID, ] %>% scale(scale=F)
# snp <- snp[rownames(snp)%in% mtb$IID, ] %>% scale(scale=F)
# if(all.equal(rownames(snp), rownames(gene)) !=T) stop('IID of snp and gene data do not match')
# mtb <- mtb[match(rownames(gene), mtb$IID), ]
# rownames(mtb) <- mtb$IID
# mtb <- mtb[,-1] %>% scale(scale=F)
# dim(gene); dim(snp); dim(mtb)
# 
# print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# r=7
# rx=10
# ry=0
# 
# print(paste(r,rx,ry))
# 
# fit_o2 <- o2m(gene, mtb, r,rx,ry)
# saveRDS(fit_o2, file=paste0('./fit_files_mtb/fit_o2_mtb_scaled_1000.rds'))
# summary(fit_o2)
# 
# fit_o2snp <- o2m(snp, mtb, r,rx,ry)
# saveRDS(fit_o2snp, file=paste0('./fit_files_mtb/fit_o2snp_mtb_scaled_1000.rds'))
# summary(fit_o2snp)
# 
# # SO2PLS with 0.01*p
# pp = 0.01
# kpx <- round(ncol(gene)*pp)
# fit_so2 <- o2m(gene, mtb, r,rx,ry, sparsity = T, keepx = kpx)
# saveRDS(fit_so2, file=paste0('./fit_files_mtb/fit_so2_mtb_scaled_1000.rds'))
# summary(fit_so2)
# 
# kpx <- round(ncol(snp)*pp)
# fit_so2snp <- o2m(snp, mtb, r,rx,ry, sparsity = T, keepx = kpx)
# saveRDS(fit_so2snp, file=paste0('./fit_files_mtb/fit_so2snp_mtb_scaled_1000.rds'))
# summary(fit_so2snp)
# 
# # GO2PLS with 0.01*groups
# SNP_gene <- readRDS('../genotypes/pgen/PCAweight_orcades.rds')
# snp_gene <- lapply(SNP_gene, function(e) data.frame(SNP=rownames(e),
#                                                     gene=rep(strsplit(colnames(e),'_')[[1]][1],nrow(e))))
# snp_gene <- bind_rows(snp_gene)
# snp_gene <- snp_gene[match(colnames(snp), snp_gene$SNP), ]
# if(all.equal(snp_gene$SNP, colnames(snp))!=T) stop('check group names')
# Nr_genes <- snp_gene$gene %>% unique %>% length
# print('Nr of genes in GO2PLS')
# print(Nr_genes)
# kpx <- round(Nr_genes*pp)
# fit_go2 <- o2m(snp, mtb, r,rx,ry, sparsity = T, groupx = snp_gene$gene, keepx = kpx, keepy = ncol(mtb))
# saveRDS(fit_go2, file=paste0('./fit_files_mtb/fit_go2_mtb_scaled_1000.rds'))
# summary(fit_go2)
# # ##################################3
# # # Scores
# # fit_o2 <- readRDS('./fit_files_mtb/fit_o2_mtb_1000.rds')
# # fit_o2snp <- readRDS('./fit_files_mtb/fit_o2snp_mtb_1000.rds')
# # fit_so2 <- readRDS('./fit_files_mtb/fit_so2_mtb_1000.rds')
# # fit_so2snp <- readRDS('./fit_files_mtb/fit_so2snp_mtb_1000.rds')
# # fit_go2 <- readRDS('./fit_files_mtb/fit_go2_mtb_1000.rds')
# # 
# # 
# T_o2_mtb <- rbind(o2_Tt(fit_o2, gene), o2_Tt(fit_o2, gene_t))
# T_o2snp_mtb <- rbind(o2_Tt(fit_o2snp, snp), o2_Tt(fit_o2snp, snp_t))
# T_so2_mtb <- rbind(o2_Tt(fit_so2, gene), o2_Tt(fit_so2, gene_t))
# T_so2snp_mtb <- rbind(o2_Tt(fit_so2snp, snp), o2_Tt(fit_so2snp, snp_t))
# T_go2_mtb <- rbind(o2_Tt(fit_go2, snp), o2_Tt(fit_go2, snp_t))
# # 
# # 
# save(T_o2_mtb, T_o2snp_mtb, T_so2_mtb, T_so2snp_mtb, T_go2_mtb, file='./fit_files_mtb/T_scaled_1000.RData')
# print('mtb unscreened Done')
# # 
# # ###########################################################################################
# # screened
# mtb <- readRDS('mtb_sel_bmi.rds')
# snp <- readRDS("../genotypes/pgen/dat_snp_orcades_mtb.rds")
# gene <- readRDS("../genotypes/pgen/GeneData_orcades_mtb.rds")
# 
# snp <- scale(snp, scale = T)
# gene <- scale(gene, scale = T)
# 
# # test data
# gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
# snp_t <- snp[rownames(snp)%in% id_test$IID, ] %>% scale(scale=F)
# dim(gene_t); dim(snp_t)
# 
# # training data
# mtb <- mtb %>% filter(IID %in% id_train$IID)
# gene <- gene[rownames(gene)%in% mtb$IID, ] %>% scale(scale=F)
# snp <- snp[rownames(snp)%in% mtb$IID, ] %>% scale(scale=F)
# if(all.equal(rownames(snp), rownames(gene)) !=T) stop('IID of snp and gene data do not match')
# mtb <- mtb[match(rownames(gene), mtb$IID), ]
# rownames(mtb) <- mtb$IID
# mtb <- mtb[,-1] %>% scale(scale=F)
# dim(gene); dim(snp); dim(mtb)
# 
# print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# # r=5
# # rx=2
# # ry=0
# r=6
# rx=5
# ry=0
# print(paste(r,rx,ry))
# 
# fit_o2 <- o2m(gene, mtb, r,rx,ry)
# saveRDS(fit_o2, file=paste0('./fit_files_mtb/fit_o2_mtb_screen_scaled_1000.rds'))
# summary(fit_o2)
# 
# fit_o2snp <- o2m(snp, mtb, r,rx,ry)
# saveRDS(fit_o2snp, file=paste0('./fit_files_mtb/fit_o2snp_mtb_screen_scaled_1000.rds'))
# summary(fit_o2snp)
# 
# # SO2PLS with 0.05*p
# pp = 0.05
# kpx <- round(ncol(gene)*pp)
# fit_so2 <- o2m(gene, mtb, r,rx,ry, sparsity = T, keepx = kpx)
# saveRDS(fit_so2, file=paste0('./fit_files_mtb/fit_so2_mtb_screen_scaled_1000.rds'))
# summary(fit_so2)
# 
# kpx <- round(ncol(snp)*pp)
# fit_so2snp <- o2m(snp, mtb, r,rx,ry, sparsity = T, keepx = kpx)
# saveRDS(fit_so2snp, file=paste0('./fit_files_mtb/fit_so2snp_mtb_screen_scaled_1000.rds'))
# summary(fit_so2snp)
# 
# # GO2PLS with 0.05*groups
# SNP_gene <- readRDS('../genotypes/pgen/PCAweight_orcades_mtb.rds')
# snp_gene <- lapply(SNP_gene, function(e) data.frame(SNP=rownames(e),
#                                                     gene=rep(strsplit(colnames(e),'_')[[1]][1],nrow(e))))
# snp_gene <- bind_rows(snp_gene)
# snp_gene <- snp_gene[match(colnames(snp), snp_gene$SNP), ]
# if(all.equal(snp_gene$SNP, colnames(snp))!=T) stop('check group names')
# Nr_genes <- snp_gene$gene %>% unique %>% length
# print('Nr of genes in GO2PLS')
# print(Nr_genes)
# kpx <- round(Nr_genes*pp)
# fit_go2 <- o2m(snp, mtb, r,rx,ry, sparsity = T, groupx = snp_gene$gene, keepx = kpx, keepy = ncol(mtb))
# saveRDS(fit_go2, file=paste0('./fit_files_mtb/fit_go2_mtb_screen_scaled_1000.rds'))
# summary(fit_go2)
# ##################################3
# # Scores
# T_o2_mtb <- rbind(o2_Tt(fit_o2, gene), o2_Tt(fit_o2, gene_t))
# T_o2snp_mtb <- rbind(o2_Tt(fit_o2snp, snp), o2_Tt(fit_o2snp, snp_t))
# T_so2_mtb <- rbind(o2_Tt(fit_so2, gene), o2_Tt(fit_so2, gene_t))
# T_so2snp_mtb <- rbind(o2_Tt(fit_so2snp, snp), o2_Tt(fit_so2snp, snp_t))
# T_go2_mtb <- rbind(o2_Tt(fit_go2, snp), o2_Tt(fit_go2, snp_t))
# 
# save(T_o2_mtb, T_o2snp_mtb, T_so2_mtb, T_so2snp_mtb, T_go2_mtb, file='./fit_files_mtb/T_screen_scaled_1000.RData')
# print('mtb screened Done')
# 

###########################################################################################
# BMI screened
mtb <- readRDS('mtb_sel_bmi.rds')
snp <- readRDS("../genotypes/pgen/dat_snp_orcades_bmi.rds")
gene <- readRDS("../genotypes/pgen/GeneData_orcades_bmi.rds")

snp <- scale(snp, scale = T)
gene <- scale(gene, scale = T)

# test data
gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
snp_t <- snp[rownames(snp)%in% id_test$IID, ] %>% scale(scale=F)
dim(gene_t); dim(snp_t)

# training data
mtb <- mtb %>% filter(IID %in% id_train$IID)
gene <- gene[rownames(gene)%in% mtb$IID, ] %>% scale(scale=F)
snp <- snp[rownames(snp)%in% mtb$IID, ] %>% scale(scale=F)
if(all.equal(rownames(snp), rownames(gene)) !=T) stop('IID of snp and gene data do not match')
mtb <- mtb[match(rownames(gene), mtb$IID), ]
rownames(mtb) <- mtb$IID
mtb <- mtb[,-1] %>% scale(scale=F)
dim(gene); dim(snp); dim(mtb)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# r=5
# rx=2
# ry=0
r=4
rx=5
ry=0
print(paste(r,rx,ry))

fit_o2 <- o2m(gene, mtb, r,rx,ry)
saveRDS(fit_o2, file=paste0('./fit_files_mtb/fit_o2_mtb_bmiscreen_scaled_1000.rds'))
summary(fit_o2)

fit_o2snp <- o2m(snp, mtb, r,rx,ry)
saveRDS(fit_o2snp, file=paste0('./fit_files_mtb/fit_o2snp_mtb_bmiscreen_scaled_1000.rds'))
summary(fit_o2snp)

# SO2PLS with 0.05*p
pp = 0.05
kpx <- round(ncol(gene)*pp)
fit_so2 <- o2m(gene, mtb, r,rx,ry, sparsity = T, keepx = kpx)
saveRDS(fit_so2, file=paste0('./fit_files_mtb/fit_so2_mtb_bmiscreen_scaled_1000.rds'))
summary(fit_so2)

kpx <- round(ncol(snp)*pp)
fit_so2snp <- o2m(snp, mtb, r,rx,ry, sparsity = T, keepx = kpx)
saveRDS(fit_so2snp, file=paste0('./fit_files_mtb/fit_so2snp_mtb_bmiscreen_scaled_1000.rds'))
summary(fit_so2snp)

# GO2PLS with 0.05*groups
SNP_gene <- readRDS('../genotypes/pgen/PCAweight_orcades_bmi.rds')
snp_gene <- lapply(SNP_gene, function(e) data.frame(SNP=rownames(e),
                                                    gene=rep(strsplit(colnames(e),'_')[[1]][1],nrow(e))))
snp_gene <- bind_rows(snp_gene)
snp_gene <- snp_gene[match(colnames(snp), snp_gene$SNP), ]
if(all.equal(snp_gene$SNP, colnames(snp))!=T) stop('check group names')
Nr_genes <- snp_gene$gene %>% unique %>% length
print('Nr of genes in GO2PLS')
print(Nr_genes)
kpx <- round(Nr_genes*pp)
fit_go2 <- o2m(snp, mtb, r,rx,ry, sparsity = T, groupx = snp_gene$gene, keepx = kpx, keepy = ncol(mtb))
saveRDS(fit_go2, file=paste0('./fit_files_mtb/fit_go2_mtb_bmiscreen_scaled_1000.rds'))
summary(fit_go2)
# ##################################3
# # Scores
# T_o2_mtb <- rbind(o2_Tt(fit_o2, gene), o2_Tt(fit_o2, gene_t))
# T_o2snp_mtb <- rbind(o2_Tt(fit_o2snp, snp), o2_Tt(fit_o2snp, snp_t))
# T_so2_mtb <- rbind(o2_Tt(fit_so2, gene), o2_Tt(fit_so2, gene_t))
# T_so2snp_mtb <- rbind(o2_Tt(fit_so2snp, snp), o2_Tt(fit_so2snp, snp_t))
# T_go2_mtb <- rbind(o2_Tt(fit_go2, snp), o2_Tt(fit_go2, snp_t))
# 
# save(T_o2_mtb, T_o2snp_mtb, T_so2_mtb, T_so2snp_mtb, T_go2_mtb, file='./fit_files_mtb/T_bmiscreen_scaled_1000.RData')
# print('mtb screened Done')

# ###########################################################################################
# # screened 5e8
# mtb <- readRDS('mtb_sel_bmi.rds')
# snp <- readRDS("../genotypes/pgen/dat_snp_orcades_mtb.rds")
# gene <- readRDS("../genotypes/pgen/GeneData_orcades_mtb5e8.rds")
# 
# snp <- scale(snp, scale = T)
# gene <- scale(gene, scale = T)
# 
# # test data
# gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
# snp_t <- snp[rownames(snp)%in% id_test$IID, ] %>% scale(scale=F)
# dim(gene_t); dim(snp_t)
# 
# # training data
# mtb <- mtb %>% filter(IID %in% id_train$IID)
# gene <- gene[rownames(gene)%in% mtb$IID, ] %>% scale(scale=F)
# snp <- snp[rownames(snp)%in% mtb$IID, ] %>% scale(scale=F)
# if(all.equal(rownames(snp), rownames(gene)) !=T) stop('IID of snp and gene data do not match')
# mtb <- mtb[match(rownames(gene), mtb$IID), ]
# rownames(mtb) <- mtb$IID
# mtb <- mtb[,-1] %>% scale(scale=F)
# dim(gene); dim(snp); dim(mtb)
# 
# print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# # r=5
# # rx=2
# # ry=0
# r=5
# rx=5
# ry=0
# print(paste(r,rx,ry))
# 
# fit_o2 <- o2m(gene, mtb, r,rx,ry)
# saveRDS(fit_o2, file=paste0('./fit_files_mtb/fit_o2_mtb_screen5e8_scaled_1000.rds'))
# summary(fit_o2)
# 
# fit_o2snp <- o2m(snp, mtb, r,rx,ry)
# saveRDS(fit_o2snp, file=paste0('./fit_files_mtb/fit_o2snp_mtb_screen5e8_scaled_1000.rds'))
# summary(fit_o2snp)
# 
# ##################################3
# # Scores
# T_o2_mtb <- rbind(o2_Tt(fit_o2, gene), o2_Tt(fit_o2, gene_t))
# T_o2snp_mtb <- rbind(o2_Tt(fit_o2snp, snp), o2_Tt(fit_o2snp, snp_t))
# 
# save(T_o2_mtb, T_o2snp_mtb, file='./fit_files_mtb/T_screen5e8_scaled_1000.RData')
# print('mtb screened Done')