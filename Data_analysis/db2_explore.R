library(dplyr)
library(moments)

id_norel <- read.table('../chr_all_QC_25.king.cutoff.in.id')

pheno <- read.table('/home/z/Data/orcades_data_for_zhujie/phenotypes/orcades_phenotypes_jay.tsv', sep = '\t', header = T)
pheno %<>% select(iid, selfrepdiab,age,sex)
pheno <- pheno %>% filter(iid %in% id_norel$V1)
pheno %<>% filter(!is.na(selfrepdiab)) %>% filter(!is.na(age)) %>% filter(!is.na(sex))
# check outlier
pheno$selfrepdiab %>% table

pheno$selfrepdiab %<>% as.factor()
pheno$sex %<>% as.factor()
pheno %>% head
colnames(pheno)[1] <- 'IID'
saveRDS(pheno, file = 'db2.rds')

id_train <- read.table('../id_train')
id_test <- read.table('../id_test')
id <- bind_rows(id_train,id_test)
(pheno %>% filter(IID %in% id$V1))$selfrepdiab %>% table

(pheno %>% filter(IID %in% id_train$V1))$selfrepdiab %>% table
