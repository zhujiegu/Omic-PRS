library(dplyr)
library(data.table)
library(OmicsPLS)
library(PO2PLS)

r=3
rx=2
ry=0
#################################
# PO2PLS with selected snps
#################################
# training

load('Intl_snp_gly_sel_training.RData')
dim(dat)
# SO2PLS
fit_snp <- PO2PLS(dat, gly, r,rx,ry, tol = 0.01, init_param = 'o2m')

saveRDS(fit_snp, file = 'fitp_snp.rds')
print(summary(fit_snp))

sd_B <- variances_inner.po2m(fit_snp, dat, gly)
z_B <- fit_snp$parameters$B %>% diag / sd_B
p_B <- 2*(1-pnorm(z_B))
# p-value of B
p_B

rm(dat); gc()
#################################
# test
load('Intl_snp_sel_test.RData')

T_sel_test <- dat_t %*% fit_snp$parameters$W
saveRDS(T_sel_test, file = 'Tp_sel_test.rds')

#################################
# O2PLS with selected snps
#################################
# # training
# 
# load('Intl_snp_gly_sel_training.RData')
# dim(dat)
# # SO2PLS
# fit_snp <- o2m(dat, gly, r,rx,ry)
# 
# saveRDS(fit_snp, file = 'fit_snp.rds')
# print(summary(fit_snp))
# rm(dat); gc()
# #################################
# # test
# load('Intl_snp_sel_test.RData')
# 
# T_sel_test <- dat_t %*% fit_snp$W.
# saveRDS(T_sel_test, file = 'T_sel_test.rds')

#################################
# SO2PLS all snps
#################################
# # training
# 
# load('Intl_snp_gly_training.RData')
# 
# # SO2PLS
# fit_snp_s <- o2m(dat, gly, r,rx,ry, sparsity = T, keepx = 100000)
# 
# saveRDS(fit_snp_s, file = 'fit_snp_s.rds')
# print(summary(fit_snp_s))
# rm(dat); gc()
# #################################
# # test
# load('Intl_snp_test.RData')
# 
# T_test <- dat_t %*% fit_snp_s$W.
# saveRDS(T_test, file = 'T_test.rds')