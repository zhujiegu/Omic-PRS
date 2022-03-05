library(dplyr)
library(data.table)
library(OmicsPLS)

r=3
rx=2
ry=0

#################################
# training

load('Intl_snp_gly_training.RData')

# SO2PLS
fit_snp_s <- o2m(dat, gly, r,rx,ry, sparsity = T, keepx = 100000)

saveRDS(fit_snp_s, file = 'fit_snp_s.rds')
print(summary(fit_snp_s))

#################################
# test