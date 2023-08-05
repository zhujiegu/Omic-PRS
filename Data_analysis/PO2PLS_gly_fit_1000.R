library(PO2PLS)
library(OmicsPLS)
library(dplyr)

id_train <- read.table('id_train')
id_test <- read.table('id_test')
colnames(id_train) <- colnames(id_test) <- 'IID'

# ######################################################################################33
print('screened scaled')

gly <- readRDS('gly_sel_bmi.rds')
gene <- readRDS("GeneData_orcades_gly.rds")

gene <- scale(gene, scale = T)

# test data
gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
dim(gene_t)

# training data
gly <- gly %>% filter(IID %in% id_train$IID)
gene <- gene[rownames(gene)%in% gly$IID, ] %>% scale(scale=F)
gly <- gly[match(rownames(gene), gly$IID), ]
rownames(gly) <- gly$IID
gly <- gly[,-1] %>% scale(scale=F)
dim(gene); dim(gly)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
r=4
rx=6
ry=0
print(paste(r,rx,ry))

fit_po2 <- PO2PLS(gene, gly, r,rx,ry, init_param = 'o2m',steps=5000, tol = 0.0001)
saveRDS(fit_po2, file=paste0('./fit_files_gly/fit_p_screen_scale_1000.rds'))

summary(fit_po2)
# fit_po2 <- readRDS('fit_po2.rds')
sd_B <- variances_inner.po2m(fit_po2, gene, gly)
fit_po2$parameters$B; sd_B
z_B <- (fit_po2$parameters$B %>% diag) / sd_B
p_B <- 2*(1-pnorm(z_B))
# p-value of B
print("p-value of global test")
print(p_B)
#########################################
# calculate scores in training and test

Tt <- with(fit_po2$parameters, (gene - gene%*%Wo%*%t(Wo)) %*% W)
Tt_t <- with(fit_po2$parameters, (gene_t - gene_t%*%Wo%*%t(Wo)) %*% W)
Tt <- rbind(Tt, Tt_t)
saveRDS(Tt, file=paste0('./fit_files_gly/T_p_screen_scale_1000.rds'))
print('Done')

######################################################################################33
# print('unscreened')
# 
# gly <- readRDS('gly_sel_bmi.rds')
# gene <- readRDS("GeneData_orcades.rds")
# 
# # test data
# gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
# dim(gene_t)
# 
# # training data
# gly <- gly %>% filter(IID %in% id_train$IID)
# gene <- gene[rownames(gene)%in% gly$IID, ] %>% scale(scale=F)
# gly <- gly[match(rownames(gene), gly$IID), ]
# rownames(gly) <- gly$IID
# gly <- gly[,-1] %>% scale(scale=F)
# dim(gene); dim(gly)
# 
# print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# r=5
# rx=2
# ry=0
# print(paste(r,rx,ry))
# 
# fit_po2 <- PO2PLS(gene, gly, r,rx,ry, init_param = 'o2m',steps=1500, tol = 0.01)
# saveRDS(fit_po2, file=paste0('./rx_gly/fit_p_',rx,'_1000.rds'))
# 
# summary(fit_po2)
# # fit_po2 <- readRDS('fit_po2.rds')
# sd_B <- variances_inner.po2m(fit_po2, gene, gly)
# fit_po2$parameters$B; sd_B
# z_B <- (fit_po2$parameters$B %>% diag) / sd_B
# p_B <- 2*(1-pnorm(z_B))
# # p-value of B
# print("p-value of global test")
# print(p_B)
# #########################################
# # calculate scores in training and test
# 
# Tt <- with(fit_po2$parameters, (gene - gene%*%Wo%*%t(Wo)) %*% W)
# Tt_t <- with(fit_po2$parameters, (gene_t - gene_t%*%Wo%*%t(Wo)) %*% W)
# Tt <- rbind(Tt, Tt_t)
# saveRDS(Tt, file=paste0('./rx_gly/T_p_',rx,'_1000.rds'))
# print('Done')

######################################################################################33
# print('BMI screened')
# 
# gly <- readRDS('gly_sel_bmi.rds')
# gene <- readRDS("../genotypes/pgen/GeneData_orcades_bmi.rds")
# gene <- scale(gene, scale = T)
# 
# # test data
# gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
# dim(gene_t)
# 
# # training data
# gly <- gly %>% filter(IID %in% id_train$IID)
# gene <- gene[rownames(gene)%in% gly$IID, ] %>% scale(scale=F)
# gly <- gly[match(rownames(gene), gly$IID), ]
# rownames(gly) <- gly$IID
# gly <- gly[,-1] %>% scale(scale=F)
# dim(gene); dim(gly)
# 
# print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# r=4
# rx=5
# ry=0
# print(paste(r,rx,ry))
# 
# fit_po2 <- PO2PLS(gene, gly, r,rx,ry, init_param = 'o2m',steps=5000, tol = 0.0001)
# saveRDS(fit_po2, file=paste0('./fit_files_gly/fit_p_bmiscreen_scaled_1000.rds'))
# 
# summary(fit_po2)
# # fit_po2 <- readRDS('fit_po2.rds')
# sd_B <- variances_inner.po2m(fit_po2, gene, gly)
# fit_po2$parameters$B; sd_B
# z_B <- (fit_po2$parameters$B %>% diag) / sd_B
# p_B <- 2*(1-pnorm(z_B))
# # p-value of B
# print("p-value of global test")
# print(p_B)
# #########################################
# # calculate scores in training and test
# 
# Tt <- with(fit_po2$parameters, (gene - gene%*%Wo%*%t(Wo)) %*% W)
# Tt_t <- with(fit_po2$parameters, (gene_t - gene_t%*%Wo%*%t(Wo)) %*% W)
# Tt <- rbind(Tt, Tt_t)
# saveRDS(Tt, file=paste0('./fit_files_gly/T_p_bmiscreen_scaled_1000.rds'))
# print('Done')