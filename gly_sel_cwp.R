library(dplyr)

gly <- readRDS('/home/z/Data/TwinsUK_CWP/Glycan.Rds')
gly <- data.frame(IID=row.names(gly), gly)

id <- read.table('id_A_decpl.txt')

cwp <- readRDS('cwp_all.rds')
colnames(cwp)[1] <- 'IID'

gly <- merge(gly, cwp, by="IID")

test <- glm(case~.-IID, data = gly, family = 'binomial')

p <- c()
# univariate
for (i in 1:76) {
  summ <- glm(paste('case~', colnames(gly)[i+1]), data = gly, family = 'binomial') %>% summary
  p[i] <- summ$coefficients[2,4]
}

p <- data.frame(IID = colnames(gly)[2:77], p=p)

gly_sel <- p %>% filter(p < 0.01)
gly_sel

##################################################################
# save glycan data with glycans associated with cwp
gly <- readRDS('/home/z/Data/TwinsUK_CWP/Glycan.Rds')
gly_sel_cwp <- gly[, colnames(gly)%in%gly_sel$IID]
saveRDS(gly_sel_cwp, file='/home/z/Data/TwinsUK_CWP/Glycan_sel.rds')

##################################################################
# change name GP to IGP
gly_sel$IID[1] <- 'IGP9'

map <- read.table('/home/z/Data/IgG_glycans_GWAS/map_r37.txt', header = T)

SNP_list <- c()
for(i in 1:nrow(gly_sel)){
  file_name <- paste('/home/z/Data/IgG_glycans_GWAS/UPLC_IgG_glycans_GWAS_', gly_sel$IID[i], '.txt.gz', sep = '')
  print(file_name)
  gwas <- read.table(gzfile(file_name), header = T)
  colnames(gwas)[1] <- 'SNP'
  gwas <- merge(gwas, map, by='SNP')
  gwas <- gwas %>% mutate(SNP=paste(Chr,':',Pos, sep = ''))
  SNP_list <- c(SNP_list, filter(gwas, p < 0.01)$SNP)
  # write.table(gwas, file=paste('/home/z/Data/IgG_glycans_GWAS/GWAS_', gly_sel$IID[i], '_edit', sep = ''), 
  #             col.names =T, quote = F, row.names = F)
}

##################################################################
# total SNP for these glycans with p-value <0.01
SNP_list %>% unique %>% length  # 93000
