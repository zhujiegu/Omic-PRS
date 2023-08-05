rx <- commandArgs(trailingOnly = T)
rx <- as.numeric(rx)

library(PO2PLS)
library(OmicsPLS)
library(dplyr)

gene <- readRDS('GeneData_training.rds')
gly <- readRDS('gly_training.rds')
if(all.equal(rownames(gly), rownames(gene))!=T){
  gly <- gly[match(rownames(gene), rownames(gly)),]
}
all.equal(rownames(gly), rownames(gene))

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
r=5
ry=0
print(paste(r,rx,ry))

fit_po2 <- PO2PLS(gene, gly, r,rx,ry, init_param = 'o2m',steps=1500, tol = 0.01)
saveRDS(fit_po2, file=paste0('./rx/fit_p_',rx,'.rds'))

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
Tt <- with(fit_po2$parameters, (gene - gene%*%Wo%*%t(Wo)) %*% W)
U <- with(fit_po2$parameters, (gly - gly%*%Co%*%t(Co)) %*% C)
jointPC <- data.frame(Tt = Tt, U = U)
To <- with(fit_po2$parameters, gene%*%Wo%*%t(Wo))

print('t % in u')
fit_po2$parameters$SigU <- with(fit_po2$parameters, SigT %*% B^2 + SigH)
with(fit_po2$parameters, SigT %*% B^2 %*% solve(SigU))

# saveRDS(jointPC, file=paste('jointPC_',r,rx,ry,'.rds', sep = ''))

cor(Tt,U)

print('Done')