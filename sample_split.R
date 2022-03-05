library(dplyr)
# library(caret)
library(magrittr)

id_A <- read.table('id_A_decpl.txt', header = T)
id_B <- read.table('id_B_decpl.txt', header = T)
id_C <- read.table('id_C_decpl.txt', header = T)


id_AB <- bind_rows(id_A, id_B)

# id_gene <- read.table('chr1.fam', header = F)
# id_gene <- read.table('chr22.psam')
# (id_AB$IID %in% id_gene$V2) %>% sum

cwp <- readRDS('cwp_all.rds')
cwp %<>% filter(PublicID %in% id_AB$IID)
colnames(cwp)[1] <- "IID"

# note the order is different
all.equal(id_AB$IID, cwp$IID)

id_AB_cwp <- merge(id_AB, cwp, by='IID')

# cwp %>% head()
id_AB_cwp %>% head
# overall proportion of cases
p <- id_AB_cwp$case %>% mean

# samplling in A
id_A_cwp <- id_AB_cwp %>% filter(IID %in% id_A$IID)
id_B_cwp <- id_AB_cwp %>% filter(IID %in% id_B$IID)

Nr_case <- as.integer(2000*p)
Nr_contl <- 2000-Nr_case

id_train <- c(sample(filter(id_A_cwp,case==1)$IID, Nr_case, replace = F),
              sample(filter(id_A_cwp,case==0)$IID, Nr_contl, replace = F))

# trainIndex <- createDataPartition(id_A_cwp$case, p = 0.706, 
#                                   list = FALSE,
#                                   times = 1)

id_A_training_cwp <- id_A_cwp %>% filter(IID %in% id_train)
id_AB_test_cwp <- bind_rows(id_A_cwp %>% filter(!IID %in% id_train), id_B_cwp)

# some check
any(id_A_training_cwp$IID %in% id_AB_test_cwp)
id_A_training_cwp$case %>% mean
id_AB_test_cwp$case %>% mean

write.table(id_A_training_cwp, file='id_A_training_cwp.txt', col.names = T, row.names = F, quote = F)
write.table(id_AB_test_cwp, file='id_AB_test_cwp.txt', col.names = T, row.names = F, quote = F)
write.table(id_AB_cwp, file='id_AB_cwp.txt', col.names = T, row.names = F, quote = F)

###################################################################
id_A_training_cwp <- read.table('id_A_training_cwp.txt', header = T)
id_AB_test_cwp <- read.table('id_AB_test_cwp.txt', header = T)

id_A_training <- id_A_training_cwp %>% select(IID)
id_AB_test <- id_AB_test_cwp %>% select(IID)
id_AB <- bind_rows(id_A_training, id_AB_test)
write.table(id_A_training, file='id_A_training.txt', col.names = T, row.names = F, quote = F)
write.table(id_AB_test, file='id_AB_test.txt', col.names = T, row.names = F, quote = F)
write.table(id_AB, file='id_AB.txt', col.names = T, row.names = F, quote = F)
