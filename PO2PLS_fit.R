library(PO2PLS)
library(OmicsPLS)
library(dplyr)

load('Intl_gene_gly_original.RData')
r=5
rx=5
ry=0

fit_po2 <- PO2PLS(gene, gly, r,rx,ry,init_param = 'o2m', tol = 0.1)
saveRDS(fit_po2, file=paste('fit_po_22gly',r,rx,ry,'.rds', sep = ''))

# fit_po2 <- readRDS('fit_po2.rds')
sd_B <- variances_inner.po2m(fit_po2, gene, gly)
z_B <- fit_po2$parameters$B %>% diag / sd_B
p_B <- 2*(1-pnorm(z_B))
# p-value of B
print("p-value of global test")
print(p_B)  # 0.000738973 0.029589419 0.199949353

#########################################
Tt <- with(fit_po2$parameters, (gene - gene%*%Wo%*%t(Wo)) %*% W)

U <- with(fit_po2$parameters, (gly - gly%*%Co%*%t(Co)) %*% C)

jointPC <- data.frame(Tt = Tt, U = U)

To <- with(fit_po2$parameters, gene%*%Wo%*%t(Wo))

print('joint % x')
ssq(Tt)/ssq(gene)
print('specific % x')
ssq(To)/ssq(gene)
print('joint % y')
ssq(U)/ssq(gly)
print('specific % y')
ssq(with(fit_po2$parameters, gly%*%Co%*%t(Co)))/ssq(gly)
print('t % in u')
fit_po2$parameters$SigU <- with(fit_po2$parameters, SigT %*% B^2 + SigH)
with(fit_po2$parameters, SigT %*% B^2 %*% solve(SigU))

# saveRDS(jointPC, file=paste('jointPC_',r,rx,ry,'.rds', sep = ''))

cor(Tt,U)
