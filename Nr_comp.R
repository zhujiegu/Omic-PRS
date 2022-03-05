library(OmicsPLS)
library(PO2PLS)

load('Intl_snp_gly_training.RData')

dim(dat)
## scree plot
fit <- o2m(dat, gly, 12, 0,0, p_thresh = 1, stripped = T)
D <- (t(fit$W.) %*% t(dat)) %*% (gly %*% fit$C.)
d <- diag(D)#[-1]
# plot(cumsum(d)/sum(d))

d1 <- svd(dat%*%t(dat), 0, 0)$d#[-1]
#plot(cumsum(d2)/sum(d2))

d2 <- svd(t(gly)%*%(gly), 0, 0)$d#[-1]
# plot(d2/sum(d2), main="Glycomics Scree plot")
# plot(cumsum(d2)/sum(d2))

save(d,d1,d2, file='eigenval_scree_snp.RData')

###############################################
# 
# load('Intl_gene_gly.RData')
# 
# ## scree plot
# fit <- o2m(gene, gly, 12, 0,0, p_thresh = 1, stripped = T)
# D <- (t(fit$W.) %*% t(gene)) %*% (gly %*% fit$C.)
# d <- diag(D)#[-1]
# # plot(cumsum(d)/sum(d))
# 
# d1 <- svd(gene%*%t(gene), 0, 0)$d#[-1]
# #plot(cumsum(d2)/sum(d2))
# 
# d2 <- svd(t(gly)%*%(gly), 0, 0)$d#[-1]
# # plot(d2/sum(d2), main="Glycomics Scree plot")
# # plot(cumsum(d2)/sum(d2))
# 
# save(d,d1,d2, file='eigenval_scree.RData')
# 
# ###############################################
load('eigenval_scree_snp.RData')

cairo_pdf(file = 'scree_plot.pdf', width = 12, height =6,
          onefile = T, fallback_resolution = 600)
par(mfrow=c(1,3))
plot(d/sum(d), main='Joint Scree plot')
plot(d1[1:100]/sum(d1[1:100]), main="Genetics Scree plot")
plot(d2/sum(d2), main="Glycomics Scree plot")
par(mfrow=c(1,1))
dev.off()
# save(d,d1,d2, file='eigenval_scree.RData')
#
# Sys.time()
# cv_nr <- vector(mode = 'list', length = 10)
# for (i in 1:10) {
#   cv_nr[[i]] <- crossval_o2m_adjR2(gene, gly, 1:6, 0:4,0:4, nr_folds = 10, nr_cores = 25)
#   print(Sys.time())
# }
#
# save(cv_nr, file='cv_nr_comp.RData')
