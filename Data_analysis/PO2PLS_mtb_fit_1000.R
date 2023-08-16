library(PO2PLS)
library(OmicsPLS)
library(dplyr)

id_train <- read.table('id_train')
id_test <- read.table('id_test')
colnames(id_train) <- colnames(id_test) <- 'IID'

######################################################################################33
print('screened scaled')

mtb <- readRDS('mtb_sel_bmi.rds')
gene <- readRDS("GeneData_orcades_mtb.rds")

gene <- scale(gene, scale = T)

# test data
gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
dim(gene_t)

# training data
mtb <- mtb %>% filter(IID %in% id_train$IID)
gene <- gene[rownames(gene)%in% mtb$IID, ] %>% scale(scale=F)
mtb <- mtb[match(rownames(gene), mtb$IID), ]
rownames(mtb) <- mtb$IID
mtb <- mtb[,-1] %>% scale(scale=F)
dim(gene); dim(mtb)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
r=6
rx=5
ry=0
print(paste(r,rx,ry))

fit_po2 <- PO2PLS(gene, mtb, r,rx,ry, init_param = 'o2m',steps=5000, tol = 0.0001)
saveRDS(fit_po2, file=paste0('./fit_files_mtb/fit_p_screen_scale_1000.rds'))

summary(fit_po2)
# fit_po2 <- readRDS('fit_po2.rds')
sd_B <- variances_inner.po2m(fit_po2, gene, mtb)
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
saveRDS(Tt, file=paste0('./fit_files_mtb/T_p_screen_scale_1000.rds'))
print('Done')
# 
# ######################################################################################33
# print('unscreened')
# 
# mtb <- readRDS('mtb_sel_bmi.rds')
# gene <- readRDS("GeneData_orcades.rds")
# 
# # test data
# gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
# dim(gene_t)
# 
# # training data
# mtb <- mtb %>% filter(IID %in% id_train$IID)
# gene <- gene[rownames(gene)%in% mtb$IID, ] %>% scale(scale=F)
# mtb <- mtb[match(rownames(gene), mtb$IID), ]
# rownames(mtb) <- mtb$IID
# mtb <- mtb[,-1] %>% scale(scale=F)
# dim(gene); dim(mtb)
# 
# print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# r=7
# rx=10
# ry=0
# print(paste(r,rx,ry))
# 
# fit_po2 <- PO2PLS(gene, mtb, r,rx,ry, init_param = 'o2m',steps=1500, tol = 0.01)
# saveRDS(fit_po2, file=paste0('./rx_mtb/fit_p_',rx,'_1000.rds'))
# 
# summary(fit_po2)
# # fit_po2 <- readRDS('fit_po2.rds')
# sd_B <- variances_inner.po2m(fit_po2, gene, mtb)
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
# saveRDS(Tt, file=paste0('./rx_mtb/T_p_',rx,'_1000.rds'))
# print('Done')

#####################################################################################33
print('BMI screened')

mtb <- readRDS('mtb_sel_bmi.rds')
gene <- readRDS("../genotypes/pgen/GeneData_orcades_bmi.rds")
gene <- scale(gene, scale = T)

# test data
gene_t <- gene[rownames(gene)%in% id_test$IID, ] %>% scale(scale=F)
dim(gene_t)

# training data
mtb <- mtb %>% filter(IID %in% id_train$IID)
gene <- gene[rownames(gene)%in% mtb$IID, ] %>% scale(scale=F)
mtb <- mtb[match(rownames(gene), mtb$IID), ]
rownames(mtb) <- mtb$IID
mtb <- mtb[,-1] %>% scale(scale=F)
dim(gene); dim(mtb)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
r=4
rx=5
ry=0
print(paste(r,rx,ry))

fit_po2 <- PO2PLS(gene, mtb, r,rx,ry, init_param = 'o2m',steps=5000, tol = 0.0001)
saveRDS(fit_po2, file=paste0('./fit_files_mtb/fit_p_bmiscreen_scaled_1000.rds'))

summary(fit_po2)
# fit_po2 <- readRDS('fit_po2.rds')
sd_B <- variances_inner.po2m(fit_po2, gene, mtb)
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
saveRDS(Tt, file=paste0('./fit_files_mtb/T_p_bmiscreen_scaled_1000.rds'))
print('Done')