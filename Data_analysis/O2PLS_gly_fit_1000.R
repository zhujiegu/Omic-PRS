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

##########################################################################################
# unscreened
gly <- readRDS('gly_sel_bmi.rds')
snp <- readRDS("../genotypes/pgen/dat_snp_orcades.rds")
gene <- readRDS("../genotypes/pgen/GeneData_orcades.rds")

snp <- scale(snp, scale = T)
gene <- scale(gene, scale = T)

# test data
gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
snp_t <- snp[rownames(snp)%in% id_test$IID, ] %>% scale(scale=F)
dim(gene_t); dim(snp_t)

# training data
gly <- gly %>% filter(IID %in% id_train$IID)
gene <- gene[rownames(gene)%in% gly$IID, ] %>% scale(scale=F)
snp <- snp[rownames(snp)%in% gly$IID, ] %>% scale(scale=F)
if(all.equal(rownames(snp), rownames(gene)) !=T) stop('IID of snp and gene data do not match')
gly <- gly[match(rownames(gene), gly$IID), ]
rownames(gly) <- gly$IID
gly <- gly[,-1] %>% scale(scale=F)
dim(gene); dim(snp); dim(gly)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
r=5
rx=2
ry=0

print(paste(r,rx,ry))

fit_o2 <- o2m(gene, gly, r,rx,ry)
saveRDS(fit_o2, file=paste0('./fit_files_gly/fit_o2_gly_scaled_1000.rds'))
summary(fit_o2)

fit_o2snp <- o2m(snp, gly, r,rx,ry)
saveRDS(fit_o2snp, file=paste0('./fit_files_gly/fit_o2snp_gly_scaled_1000.rds'))
summary(fit_o2snp)

# SO2PLS with 0.01*p
pp = 0.001
kpx <- round(ncol(gene)*pp)
fit_so2 <- o2m(gene, gly, r,rx,ry, sparsity = T, keepx = kpx)
saveRDS(fit_so2, file=paste0('./fit_files_gly/fit_so2_gly_scaled_1000.rds'))
summary(fit_so2)

kpx <- round(ncol(snp)*pp)
fit_so2snp <- o2m(snp, gly, r,rx,ry, sparsity = T, keepx = kpx)
saveRDS(fit_so2snp, file=paste0('./fit_files_gly/fit_so2snp_gly_scaled_1000.rds'))
summary(fit_so2snp)

# GO2PLS with 0.01*groups
SNP_gene <- readRDS('../genotypes/pgen/PCAweight_orcades.rds')
snp_gene <- lapply(SNP_gene, function(e) data.frame(SNP=rownames(e),
                                                    gene=rep(strsplit(colnames(e),'_')[[1]][1],nrow(e))))
snp_gene <- bind_rows(snp_gene)
snp_gene <- snp_gene[match(colnames(snp), snp_gene$SNP), ]
if(all.equal(snp_gene$SNP, colnames(snp))!=T) stop('check group names')
Nr_genes <- snp_gene$gene %>% unique %>% length
print('Nr of genes in GO2PLS')
print(Nr_genes)
kpx <- round(Nr_genes*pp)
fit_go2 <- o2m(snp, gly, r,rx,ry, sparsity = T, groupx = snp_gene$gene, keepx = kpx, keepy = ncol(gly))
saveRDS(fit_go2, file=paste0('./fit_files_gly/fit_go2_gly_scaled_1000.rds'))
summary(fit_go2)
#
# ##################################3
# # Scores
# fit_o2 <- readRDS('./fit_files_gly/fit_o2_gly_1000.rds')
# fit_o2snp <- readRDS('./fit_files_gly/fit_o2snp_gly_1000.rds')
# fit_so2 <- readRDS('./fit_files_gly/fit_so2_gly_1000.rds')
# fit_so2snp <- readRDS('./fit_files_gly/fit_so2snp_gly_1000.rds')
# fit_go2 <- readRDS('./fit_files_gly/fit_go2_gly_1000.rds')
#
#
T_o2_gly <- rbind(o2_Tt(fit_o2, gene), o2_Tt(fit_o2, gene_t))
T_o2snp_gly <- rbind(o2_Tt(fit_o2snp, snp), o2_Tt(fit_o2snp, snp_t))
T_so2_gly <- rbind(o2_Tt(fit_so2, gene), o2_Tt(fit_so2, gene_t))
T_so2snp_gly <- rbind(o2_Tt(fit_so2snp, snp), o2_Tt(fit_so2snp, snp_t))
T_go2_gly <- rbind(o2_Tt(fit_go2, snp), o2_Tt(fit_go2, snp_t))

save(T_o2_gly, T_o2snp_gly, T_so2_gly, T_so2snp_gly, T_go2_gly, file='./fit_files_gly/T_scaled_1000.RData')
print('gly unscreened Done')

###########################################################################################
# screened
gly <- readRDS('gly_sel_bmi.rds')
snp <- readRDS("../genotypes/pgen/dat_snp_orcades_gly.rds")
gene <- readRDS("../genotypes/pgen/GeneData_orcades_gly.rds")

snp <- scale(snp, scale = T)
gene <- scale(gene, scale = T)
# test data
gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
snp_t <- snp[rownames(snp)%in% id_test$IID, ] %>% scale(scale=F)
dim(gene_t); dim(snp_t)

# training data
gly <- gly %>% filter(IID %in% id_train$IID)
gene <- gene[rownames(gene)%in% gly$IID, ] %>% scale(scale=F)
snp <- snp[rownames(snp)%in% gly$IID, ] %>% scale(scale=F)
if(all.equal(rownames(snp), rownames(gene)) !=T) stop('IID of snp and gene data do not match')
gly <- gly[match(rownames(gene), gly$IID), ]
rownames(gly) <- gly$IID
gly <- gly[,-1] %>% scale(scale=F)
dim(gene); dim(snp); dim(gly)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# r=5
# rx=2
# ry=0
r=4
rx=3
ry=0
print(paste(r,rx,ry))

fit_o2 <- o2m(gene, gly, r,rx,ry)
saveRDS(fit_o2, file=paste0('./fit_files_gly/fit_o2_gly_screen_scaled_1000.rds'))
summary(fit_o2)

fit_o2snp <- o2m(snp, gly, r,rx,ry)
saveRDS(fit_o2snp, file=paste0('./fit_files_gly/fit_o2snp_gly_screen_scaled_1000.rds'))
summary(fit_o2snp)

# SO2PLS with 0.05*p
pp = 0.05
kpx <- round(ncol(gene)*pp)
fit_so2 <- o2m(gene, gly, r,rx,ry, sparsity = T, keepx = kpx)
saveRDS(fit_so2, file=paste0('./fit_files_gly/fit_so2_gly_screen_scaled_1000.rds'))
summary(fit_so2)

kpx <- round(ncol(snp)*pp)
fit_so2snp <- o2m(snp, gly, r,rx,ry, sparsity = T, keepx = kpx)
saveRDS(fit_so2snp, file=paste0('./fit_files_gly/fit_so2snp_gly_screen_scaled_1000.rds'))
summary(fit_so2snp)

# GO2PLS with 0.05*groups
SNP_gene <- readRDS('../genotypes/pgen/PCAweight_orcades_gly.rds')
snp_gene <- lapply(SNP_gene, function(e) data.frame(SNP=rownames(e),
                                                    gene=rep(strsplit(colnames(e),'_')[[1]][1],nrow(e))))
snp_gene <- bind_rows(snp_gene)
snp_gene <- snp_gene[match(colnames(snp), snp_gene$SNP), ]
if(all.equal(snp_gene$SNP, colnames(snp))!=T) stop('check group names')
Nr_genes <- snp_gene$gene %>% unique %>% length
print('Nr of genes in GO2PLS')
print(Nr_genes)
kpx <- round(Nr_genes*pp)
fit_go2 <- o2m(snp, gly, r,rx,ry, sparsity = T, groupx = snp_gene$gene, keepx = kpx, keepy = ncol(gly))
saveRDS(fit_go2, file=paste0('./fit_files_gly/fit_go2_gly_screen_scaled_1000.rds'))
summary(fit_go2)

##################################3
# Scores
T_o2_gly <- rbind(o2_Tt(fit_o2, gene), o2_Tt(fit_o2, gene_t))
T_o2snp_gly <- rbind(o2_Tt(fit_o2snp, snp), o2_Tt(fit_o2snp, snp_t))
T_so2_gly <- rbind(o2_Tt(fit_so2, gene), o2_Tt(fit_so2, gene_t))
T_so2snp_gly <- rbind(o2_Tt(fit_so2snp, snp), o2_Tt(fit_so2snp, snp_t))
T_go2_gly <- rbind(o2_Tt(fit_go2, snp), o2_Tt(fit_go2, snp_t))

save(T_o2_gly, T_o2snp_gly, T_so2_gly, T_so2snp_gly, T_go2_gly, file='./fit_files_gly/T_screen_scaled_1000.RData')
print('gly screened Done')



###########################################################################################
# BMI screened
gly <- readRDS('gly_sel_bmi.rds')
snp <- readRDS("../genotypes/pgen/dat_snp_orcades_bmi.rds")
gene <- readRDS("../genotypes/pgen/GeneData_orcades_bmi.rds")

snp <- scale(snp, scale = T)
gene <- scale(gene, scale = T)

# test data
gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
snp_t <- snp[rownames(snp)%in% id_test$IID, ] %>% scale(scale=F)
dim(gene_t); dim(snp_t)

# training data
gly <- gly %>% filter(IID %in% id_train$IID)
gene <- gene[rownames(gene)%in% gly$IID, ] %>% scale(scale=F)
snp <- snp[rownames(snp)%in% gly$IID, ] %>% scale(scale=F)
if(all.equal(rownames(snp), rownames(gene)) !=T) stop('IID of snp and gene data do not match')
gly <- gly[match(rownames(gene), gly$IID), ]
rownames(gly) <- gly$IID
gly <- gly[,-1] %>% scale(scale=F)
dim(gene); dim(snp); dim(gly)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
r=4
rx=5
ry=0
print(paste(r,rx,ry))

fit_o2 <- o2m(gene, gly, r,rx,ry)
saveRDS(fit_o2, file=paste0('./fit_files_gly/fit_o2_gly_bmiscreen_scaled_1000.rds'))
summary(fit_o2)

fit_o2snp <- o2m(snp, gly, r,rx,ry)
saveRDS(fit_o2snp, file=paste0('./fit_files_gly/fit_o2snp_gly_bmiscreen_scaled_1000.rds'))
summary(fit_o2snp)

# SO2PLS with 0.05*p
pp = 0.05
kpx <- round(ncol(gene)*pp)
fit_so2 <- o2m(gene, gly, r,rx,ry, sparsity = T, keepx = kpx)
saveRDS(fit_so2, file=paste0('./fit_files_gly/fit_so2_gly_bmiscreen_scaled_1000.rds'))
summary(fit_so2)

kpx <- round(ncol(snp)*pp)
fit_so2snp <- o2m(snp, gly, r,rx,ry, sparsity = T, keepx = kpx)
saveRDS(fit_so2snp, file=paste0('./fit_files_gly/fit_so2snp_gly_bmiscreen_scaled_1000.rds'))
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
fit_go2 <- o2m(snp, gly, r,rx,ry, sparsity = T, groupx = snp_gene$gene, keepx = kpx, keepy = ncol(gly))
saveRDS(fit_go2, file=paste0('./fit_files_gly/fit_go2_gly_bmiscreen_scaled_1000.rds'))
summary(fit_go2)

##################################3
# Scores
T_o2_gly <- rbind(o2_Tt(fit_o2, gene), o2_Tt(fit_o2, gene_t))
T_o2snp_gly <- rbind(o2_Tt(fit_o2snp, snp), o2_Tt(fit_o2snp, snp_t))
T_so2_gly <- rbind(o2_Tt(fit_so2, gene), o2_Tt(fit_so2, gene_t))
T_so2snp_gly <- rbind(o2_Tt(fit_so2snp, snp), o2_Tt(fit_so2snp, snp_t))
T_go2_gly <- rbind(o2_Tt(fit_go2, snp), o2_Tt(fit_go2, snp_t))

save(T_o2_gly, T_o2snp_gly, T_so2_gly, T_so2snp_gly, T_go2_gly, file='./fit_files_gly/T_bmiscreen_scaled_1000.RData')
print('gly screened Done')


#########################################################################################
# screened 5e-8
gly <- readRDS('gly_sel_bmi.rds')
snp <- readRDS("../genotypes/pgen/dat_snp_orcades_gly.rds")
gene <- readRDS("../genotypes/pgen/GeneData_orcades_gly5e8.rds")

snp <- scale(snp, scale = T)
gene <- scale(gene, scale = T)
# test data
gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
snp_t <- snp[rownames(snp)%in% id_test$IID, ] %>% scale(scale=F)
dim(gene_t); dim(snp_t)

# training data
gly <- gly %>% filter(IID %in% id_train$IID)
gene <- gene[rownames(gene)%in% gly$IID, ] %>% scale(scale=F)
snp <- snp[rownames(snp)%in% gly$IID, ] %>% scale(scale=F)
if(all.equal(rownames(snp), rownames(gene)) !=T) stop('IID of snp and gene data do not match')
gly <- gly[match(rownames(gene), gly$IID), ]
rownames(gly) <- gly$IID
gly <- gly[,-1] %>% scale(scale=F)
dim(gene); dim(snp); dim(gly)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# r=5
# rx=2
# ry=0
r=4
rx=5
ry=0
print(paste(r,rx,ry))

fit_o2 <- o2m(gene, gly, r,rx,ry)
saveRDS(fit_o2, file=paste0('./fit_files_gly/fit_o2_gly_screen5e8_scaled_1000.rds'))
summary(fit_o2)

fit_o2snp <- o2m(snp, gly, r,rx,ry)
saveRDS(fit_o2snp, file=paste0('./fit_files_gly/fit_o2snp_gly_screen5e8_scaled_1000.rds'))
summary(fit_o2snp)

##################################3
# Scores
T_o2_gly <- rbind(o2_Tt(fit_o2, gene), o2_Tt(fit_o2, gene_t))
T_o2snp_gly <- rbind(o2_Tt(fit_o2snp, snp), o2_Tt(fit_o2snp, snp_t))

save(T_o2_gly, T_o2snp_gly, file='./fit_files_gly/T_screen5e8_scaled_1000.RData')
print('gly screened Done')