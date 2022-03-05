library(magrittr)
library(dplyr)

psam <- read.table('chr1.psam')
colnames(psam)[1] <- 'IID'
id_A_training <- read.table('id_A_training.txt', header = T)

# header + number of lines
row_training <- c(1, which(psam$IID %in% id_A_training$IID)+1)
row_test <- c(1, which(!psam$IID %in% id_A_training$IID)+1)

write.table(row_training, file='row_training', col.names = F, row.names = F, quote = F)
write.table(row_test, file='row_test', col.names = F, row.names = F, quote = F)

psam_training <- psam %>% filter(IID %in% id_A_training$IID) %>% select(IID)
psam_test <- psam %>% filter(!IID %in% id_A_training$IID) %>% select(IID)

write.table(psam_training, file='psam_training', col.names = T, row.names = F, quote = F)
write.table(psam_test, file='psam_test', col.names = T, row.names = F, quote = F)
