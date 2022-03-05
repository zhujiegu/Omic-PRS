library(dplyr)
library(magrittr)

source('/home/z/Mygit/GLM_PO2PLS/glmPO2PLS.R')

load('Intl_gene_gly.RData')
gly <- readRDS('/home/z/Data/TwinsUK_CWP/Glycan.Rds')

jointPC <- readRDS('jointPC_490.rds')
fit_o2 <- readRDS('fit_o2_490.rds')
fit_po <- readRDS('fit_po490.rds')

cor(fit_po$parameters$W, fit_o2$W.)
cor(fit_po$parameters$Wo, fit_o2$W.)
cor(fit_po$parameters$Wo, fit_o2$W_Yosc)

cor(jointPC, fit_o2$Tt)
cor(jointPC, fit_o2$U)

cor(jointPC, fit_o2$T_Yosc)

# PCA
gen_pca <- svd(gene, nu = 0, nv = 20)

# T_PCA <- gene%*%gen_pca$v
# cor(T_PCA, jointPC %>% dplyr::select(starts_with('T')))
cor(gen_pca$v, fit_po$parameters$Wo) %>% round(digits = 2)
cor(gen_pca$v, fit_po$parameters$W) %>% round(digits = 2)

gly_pca <- svd(gly, nu = 0, nv = 20)
cor(gly_pca$v, fit_po$parameters$C) %>% round(digits = 2)

cor(gly_pca$v, fit_o2$C.) %>% round(digits = 2)
