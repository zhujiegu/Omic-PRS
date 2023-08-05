library(dplyr)
library(data.table)

map <- fread('/home/z/Data/IgG_glycans_GWAS/map_r37.txt', header = T)

# SNP_list <- c()
for(i in 1:23){
  file_name <- paste0('/home/z/Data/IgG_glycans_GWAS/UPLC_IgG_glycans_GWAS_IGP', i, '.txt.gz')
  print(file_name)
  gwas <- read.table(gzfile(file_name), header = T)
  colnames(gwas)[1] <- 'SNP'
  gwas <- merge(gwas, map, by='SNP')
  gwas <- gwas %>% mutate(SNPpos=paste(Chr,':',Pos, sep = ''))
  # SNP_list <- c(SNP_list, filter(gwas, p < 0.01)$SNP)
  fwrite(gwas, file=paste0('/home/z/Data/IgG_glycans_GWAS/GWAS_IGP', i, '_edit'),
              col.names =T, quote = F, row.names = F, sep = '\t')
}

##################################################################
# # total SNP for these glycans with p-value <0.01
# SNP_list <- SNP_list %>% unique
# SNP_list %>% length  # 93000
# saveRDS(SNP_list, file = 'SNPs_gly.rds')
