library(dplyr)
library(magrittr)

gly <- read.table(file = '/home/z/Data/orcades_data_for_zhujie/glycomics/orcades_IgG_QC_data_NEW.tsv', sep = '\t', header = TRUE)
pheno <- read.table('/home/z/Data/orcades_data_for_zhujie/phenotypes/orcades_phenotypes_jay.tsv', sep = '\t', header = T)

# original peaks
gly <- gly[,1:24]
gly %>% head

gly[,-1] <- log10(gly[,-1]) %>% scale(scale = F)
gly[,2:24] %>% boxplot

# correct for age and sex
covar <- pheno %>% select(iid, age, sex)
covar$sex %<>% as.factor
gly <- merge(gly, covar, by='iid')
gly$age %>% is.na %>% sum
gly$sex %>% is.na %>% sum

gly$age %>% density %>% plot
gly$sex %>% table

correct <- lapply(colnames(gly)[2:24], function(j){
  fit <- lm(paste0(j,'~age+sex'), data = gly)
  return(residuals(fit))
})

correct_matrix <- do.call(cbind, correct)
gly[,2:24] <- correct_matrix
gly <- gly %>% select(-age,-sex)

saveRDS(gly, file='gly_gp_orcades.rds')
