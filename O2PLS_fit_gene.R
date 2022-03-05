library(OmicsPLS)
library(dplyr)

load('Intl_gene_gly_sel.RData')
r=3
rx=9
ry=0

#############################
fit_s <- o2m(gene, gly, r,rx,ry, sparsity = T, keepx = 1000)

print(summary(fit_s))

glm(id$case~fit_s$Tt, family = 'binomial') %>% summary
saveRDS(fit_s, file = 'fit_s_100.rds')

#############################

load('Gene_cwp_test.RData')
fit_s_100 <- readRDS('fit_s_100.rds')

Tt_test <- gene_test %*% fit_s_100$W.
glm(id_test$case~Tt_test, family = 'binomial') %>% summary


fit_o2 <- readRDS('fit_o2.rds')
cor(fit_o2$Tt, fit_o2$U)

summary(fit_o2)

cor(Tt, fit_o2$Tt)
