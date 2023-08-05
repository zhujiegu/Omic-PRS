x <- readRDS('dat_snp.rds')
gene <- readRDS('GeneData.rds')

pca <- prcomp(gene, rank.=20)
saveRDS(pca, file = 'PCA_gene.rds')

pca_snp <- prcomp(x, rank.=20)
saveRDS(pca_snp, file = 'PCA_snp.rds')

print('PCA of x done')