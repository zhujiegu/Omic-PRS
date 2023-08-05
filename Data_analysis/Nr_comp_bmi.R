library(OmicsPLS)
library(PO2PLS)
library(dplyr)


gly <- readRDS('gly_sel_bmi.rds')
gene <- readRDS("../genotypes/pgen/GeneData_orcades_bmi.rds")

gene <- gene[rownames(gene)%in% gly$IID, ] %>% scale(scale=F)
gly <- gly[match(rownames(gene), gly$IID), ]
rownames(gly) <- gly$IID
gly <- gly[,-1] %>% scale(scale=F)
all.equal(rownames(gene), rownames(gly))

## scree plot
fit <- o2m(gene, gly, 10, 0,0, p_thresh = 1, stripped = T)
D <- (t(fit$W.) %*% t(gene)) %*% (gly %*% fit$C.)
d <- diag(D)#[-1]
# plot(cumsum(d)/sum(d))

d1 <- svd(gene%*%t(gene), 0, 0)$d#[-1]
#plot(cumsum(d2)/sum(d2))

d2 <- svd(t(gly)%*%(gly), 0, 0)$d#[-1]
# plot(d2/sum(d2), main="Glycomics Scree plot")
# plot(cumsum(d2)/sum(d2))

save(d,d1,d2, file='eigenval_scree_gly_bmiscreen.RData')

#####################################################################
# load('eigenval_scree_gly.RData')

cairo_pdf(file = 'scree_plot_gly_bmiscreen.pdf', width = 12, height =6,
          onefile = T, fallback_resolution = 600)
par(mfrow=c(1,3))
plot(d/sum(d), main='Joint Scree plot')
plot(d1[1:100]/sum(d1[1:100]), main="Genetics Scree plot")
plot(d2/sum(d2), main="Glycomics Scree plot")
par(mfrow=c(1,1))
dev.off()
# cumsum(d1)/sum(d1) > 0.1
#####################################################################
gly <- readRDS('gly_sel_bmi.rds')
snp <- readRDS("../genotypes/pgen/dat_snp_orcades_bmi.rds")

snp <- snp[rownames(snp)%in% gly$IID, ] %>% scale(scale=F)
gly <- gly[match(rownames(snp), gly$IID), ]
rownames(gly) <- gly$IID
gly <- gly[,-1] %>% scale(scale=F)
all.equal(rownames(snp), rownames(gly))

## scree plot
fit <- o2m(snp, gly, 10, 0,0, p_thresh = 1, stripped = T)
D <- (t(fit$W.) %*% t(snp)) %*% (gly %*% fit$C.)
d <- diag(D)#[-1]
# plot(cumsum(d)/sum(d))

d1 <- svd(snp%*%t(snp), 0, 0)$d#[-1]
#plot(cumsum(d2)/sum(d2))

d2 <- svd(t(gly)%*%(gly), 0, 0)$d#[-1]
# plot(d2/sum(d2), main="Glycomics Scree plot")
# plot(cumsum(d2)/sum(d2))

save(d,d1,d2, file='eigenval_scree_snp_gly_bmiscreen.RData')

###############################
# load('eigenval_scree_snp_gly.RData')

cairo_pdf(file = 'scree_plot_snp_gly_bmiscreen.pdf', width = 12, height =6,
          onefile = T, fallback_resolution = 600)
par(mfrow=c(1,3))
plot(d/sum(d), main='Joint Scree plot')
plot(d1[1:100]/sum(d1[1:100]), main="Genetics Scree plot")
plot(d2/sum(d2), main="Glycomics Scree plot")
par(mfrow=c(1,1))
dev.off()
# cumsum(d1)/sum(d1) > 0.1

#############################################################################################
mtb <- readRDS('mtb_sel_bmi.rds')
gene <- readRDS("../genotypes/pgen/GeneData_orcades_bmi.rds")

gene <- gene[rownames(gene)%in% mtb$IID, ] %>% scale(scale=F)
mtb <- mtb[match(rownames(gene), mtb$IID), ]
rownames(mtb) <- mtb$IID
mtb <- mtb[,-1] %>% scale(scale=F)
all.equal(rownames(gene), rownames(mtb))

## scree plot
fit <- o2m(gene, mtb, 50, 0,0, p_thresh = 1, stripped = T)
D <- (t(fit$W.) %*% t(gene)) %*% (mtb %*% fit$C.)
d <- diag(D)#[-1]
# plot(cumsum(d)/sum(d))

d1 <- svd(gene%*%t(gene), 0, 0)$d#[-1]
#plot(cumsum(d2)/sum(d2))

d2 <- svd(t(mtb)%*%(mtb), 0, 0)$d#[-1]
# plot(d2/sum(d2), main="mtbcomics Scree plot")
# plot(cumsum(d2)/sum(d2))

save(d,d1,d2, file='eigenval_scree_mtb_bmiscreen.RData')
print('Gene done')
#####################################################################
# load('eigenval_scree_mtb.RData')

cairo_pdf(file = 'scree_plot_mtb_bmiscreen.pdf', width = 12, height =6,
          onefile = T, fallback_resolution = 600)
par(mfrow=c(1,3))
plot(d/sum(d), main='Joint Scree plot')
plot(d1[1:100]/sum(d1[1:100]), main="Genetics Scree plot")
plot(d2/sum(d2), main="mtbcomics Scree plot")
par(mfrow=c(1,1))
dev.off()
# cumsum(d1)/sum(d1) > 0.1

#####################################################################
mtb <- readRDS('mtb_sel_bmi.rds')
snp <- readRDS("../genotypes/pgen/dat_snp_orcades_bmi.rds")

snp <- snp[rownames(snp)%in% mtb$IID, ] %>% scale(scale=F)
mtb <- mtb[match(rownames(snp), mtb$IID), ]
rownames(mtb) <- mtb$IID
mtb <- mtb[,-1] %>% scale(scale=F)
all.equal(rownames(snp), rownames(mtb))

## scree plot
fit <- o2m(snp, mtb, 50, 0,0, p_thresh = 1, stripped = T)
D <- (t(fit$W.) %*% t(snp)) %*% (mtb %*% fit$C.)
d <- diag(D)#[-1]
# plot(cumsum(d)/sum(d))

d1 <- svd(snp%*%t(snp), 0, 0)$d#[-1]
#plot(cumsum(d2)/sum(d2))

d2 <- svd(t(mtb)%*%(mtb), 0, 0)$d#[-1]
# plot(d2/sum(d2), main="mtbcomics Scree plot")
# plot(cumsum(d2)/sum(d2))

save(d,d1,d2, file='eigenval_scree_snp_mtb_bmiscreen.RData')
print('snp done')


###############################
# load('eigenval_scree_snp_mtb.RData')

cairo_pdf(file = 'scree_plot_snp_mtb_bmiscreen.pdf', width = 12, height =6,
          onefile = T, fallback_resolution = 600)
par(mfrow=c(1,3))
plot(d/sum(d), main='Joint Scree plot')
plot(d1[1:100]/sum(d1[1:100]), main="Genetics Scree plot")
plot(d2/sum(d2), main="mtbcomics Scree plot")
par(mfrow=c(1,1))
dev.off()
# cumsum(d1)/sum(d1) > 0.1

