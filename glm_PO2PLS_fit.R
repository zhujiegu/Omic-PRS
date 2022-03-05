library(dplyr)
library(magrittr)
library(PO2PLS)
library(OmicsPLS)

source('glmPO2PLS.R')
# source('/home/z/Mygit/GLM_PO2PLS/glmPO2PLS.R')

load('Intl_gene_gly.RData')
cwp <- scale(cwp, scale = F)

r=6
rx=7
ry=3

fit <- glm_PO2PLS(gene, gly, cwp, r,rx,ry, init_param = 'o2m', tol = 0.1, family = 'Gaussian')
saveRDS(fit, file=paste('fit_glmp',r,rx,ry,'.rds', sep = ''))

sd_B <- sd_B(fit, gene, gly, cwp)
z_B <- fit$params$B %>% diag / sd_B
p_B <- 2*(1-pnorm(z_B))
# p-value of B
print("p-value of global test")
print(p_B)  # 0.000738973 0.029589419 0.199949353

#########################################
Tt <- with(fit$params, (gene - gene%*%Wo%*%t(Wo)) %*% W)

U <- with(fit$params, (gly - gly%*%Co%*%t(Co)) %*% C)

jointPC <- data.frame(Tt = Tt, U = U)

To <- with(fit$params, gene%*%Wo%*%t(Wo))

print('joint % x')
ssq(Tt)/ssq(gene)
print('specific % x')
ssq(To)/ssq(gene)
print('joint % y')
ssq(U)/ssq(gly)
print('specific % y')
ssq(with(fit$params, gly%*%Co%*%t(Co)))/ssq(gly)
print('t % in u')
with(fit$params, SigT %*% B^2 %*% solve(SigU))

saveRDS(jointPC, file=paste('jointPC_',r,rx,ry,'.rds', sep = ''))

cor(Tt,U)
