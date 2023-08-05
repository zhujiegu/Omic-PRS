library(dplyr)

id_gene <- read.table('chr_all_QC_25.king.cutoff.in.id')
id_gene <- id_gene[,1,drop=F]

bmi <- readRDS('./BMI/bmi.rds')
gly <- readRDS('./BMI/gly_sel_bmi.rds')
mtb <- readRDS('./BMI/mtb_sel_bmi.rds')

id_gene %>% filter(V1 %in% gly$IID) %>% filter(V1 %in% mtb$IID) %>% filter(!V1 %in% bmi$IID) %>% nrow

# overlapping all
id_all <- id_gene %>% filter(V1 %in% gly$IID) %>% filter(V1 %in% mtb$IID) %>% 
  filter(V1 %in% bmi$IID) %>% select(V1)
colnames(id_all) <- 'IID'

train_idx <- sample(1:nrow(id_all), 1000, replace = F)
id_train <- id_all[train_idx,]
id_test <- id_all[-train_idx,]

write.table(id_train, file = 'id_train', col.names = F, row.names = F, quote = F)
write.table(id_test, file = 'id_test', col.names = F, row.names = F, quote = F)

# BD II
db <- readRDS('db2.rds')
db %>% filter(IID %in% id_train$IID) %>% select(selfrepdiab) %>% table
db %>% filter(IID %in% id_test$IID) %>% select(selfrepdiab) %>% table
