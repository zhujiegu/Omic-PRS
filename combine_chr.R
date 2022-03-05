library(dplyr)
library(data.table)

#####################################################
# SNP data .raw files of plink output
chr <- vector(mode = "list", length = 22)

# too slow, use cut and paste directly
# for(i in 1:22){
#   print(i)
#   # chr[[i]] <- fread(paste('chr',i,'.raw', sep = ''), header = T, verbose = F, select = 10:20)
#   chr[[i]] <- fread(paste('chr',i,'.raw', sep = ''), header = T, drop = 1:6, verbose = F)
# }
# 
# dat <- do.call("cbind", chr)
# rm(chr); gc()
sam_id <- read.table('chr1.psam')
rownames(dat) <- sam_id$V1

snp_id <- colnames(dat)

snp <- sapply(snp_id, function(e) strsplit(e,'_')[[1]][1])
colnames(dat) <- snp
print(dim(dat))

print(dat[1:6,1:6])

id_train <- read.table('id_A_training_cwp.txt', header = T)
id_test <- read.table('id_AB_test_cwp.txt', header = T)

dat_train <- dat[match(id_train$IID, rownames(dat)), ]
dat_test <- dat[match(id_test$IID, rownames(dat)), ]

all.equal(rownames(dat_train), id_train$IID) %>% print
all.equal(rownames(dat_test), id_test$IID) %>% print

saveRDS(dat_train, file='SNP_dat_train.rds')
saveRDS(dat_test, file='SNP_dat_test.rds')
#####################################################
# # summarized gene data
# 
# # load data of each chr
# files <- list.files(pattern = "Gene.+.rds")
# dat <- do.call("cbind", lapply(files, readRDS))
# 
# # load sample id
# sam_id <- read.table('chr1.psam')
# 
# dim(dat)
# length(sam_id$V1)
# 
# rownames(dat) <- sam_id$V1
# 
# sapply(strsplit(colnames(dat), "_"), function(e) e[1]) %>% unique %>% length
# 
# saveRDS(dat, 'GeneDat.Rds')


