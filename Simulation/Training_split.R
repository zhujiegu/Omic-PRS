library(dplyr)
# library(caret)
library(magrittr)

arrID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(arrID)

id <- read.table('id_AB.txt', header = T)

id_train <- sample(id$IID, 2000, replace = F)
id_train <- data.frame(IID=id_train)

write.table(id_train, file=paste('./dat_rep/id_training_',arrID,'.txt', sep = ''), col.names = T, row.names = F, quote = F)

###################
# # Nov 1 /2022 create tune and target in test
# id <- read.table('id_AB.txt', header = T)
# for (arrID in 1:50) {
#   id_train <- read.table(file=paste0('./dat_rep/id_training_',arrID,'.txt'), header = T)
#   id_test <- id %>% filter(!IID %in% id_train$IID)
#   set.seed(arrID)
#   id_tune <- sample(id_test$IID, 465, replace = F)
#   id_tune <- data.frame(IID=id_tune)
#   id_target <- id_test %>% filter(!IID %in% id_tune$IID)
#   write.table(id_tune, file=paste('./dat_rep/id_tune_',arrID,'.txt', sep = ''), 
#               col.names = T, row.names = F, quote = F)
#   write.table(id_target, file=paste('./dat_rep/id_target_',arrID,'.txt', sep = ''), 
#               col.names = T, row.names = F, quote = F)
# }

