library(data.table)
library(lassosum)
library(dplyr)
library(parallel)

id_snp <- read.table('chr1.fam')
id_snp <- id_snp$V2

id_A <- read.table('id_A_decpl.txt', header = T)
id_B <- read.table('id_B_decpl.txt', header = T)
id_C <- read.table('id_C_decpl.txt', header = T)

# reference panal from (C +1500 from A)
draw_ref <- sample(id_A$IID, 1500, replace = F)
id_ref <- c(id_C$IID, draw_ref)

# lassosum ref.keep format
keep_ref <- id_snp %in% id_ref

### Read summary statistics file ###
ss <- fread("../CWP_GWAS_EU_ANCESTRY_UKB.txt")
head(ss)

### get cor of X and Y from p-value
cor <- p2cor(p = ss$P, n = 249843)

### Read LD region file ###
LDblocks <- "EUR.hg19"

N=length(id_snp)

out <- vector(mode = "list", length = 22)
for (k in 1:22) {
  # here, the 1000 sample are used as reference panel, 
  # 500 as test (or, actually validation data for tuning hyperparameters
  # rest as the target, or the independent test
  # 8T of memory needed

  ### Specify the PLINK file stub of the reference panel ###
  ref.bfile <- paste('chr',k, sep = '')
  
  ### Specify the PLINK file stub of the test data ###
  test.bfile <- paste('chr',k, sep = '')
  
  out[[k]] <- lassosum.pipeline(cor=cor, chr=ss$CHR, pos=ss$BP, 
                           A1=ss$A1, A2=ss$A2, # A2 is not required but advised
                           ref.bfile=ref.bfile, test.bfile=test.bfile, 
                           LDblocks = LDblocks, keep.ref = keep_ref, 
                           remove.test = keep_ref)
}

#
# out <- readRDS('out_lassosum.rds')
out_all <- readRDS('out_all_lassosum.rds')
# combine chrs
out_all <- do.call(merge, out)

# 500 samples to select hyperparameters
id_val <- id_snp[!id_snp %in% id_ref] # not in ref
id_val <- id_val[id_val %in% id_A$IID] # in Block B
draw_val <- sample(1:length(id_val), 500, replace = F) %>% sort() # consider improve by keeping case/control ratio
id_val <- id_val[draw_val]
keep_val <- id_snp %in% id_val
all.equal(id_val, id_snp[keep_val])

### outcome and covariates
cwp <- readRDS('cwp_all.rds')
covar <- readRDS('Covar.Rds')
# cwp <- cwp %>% filter(cwp$PublicID %in% id_snp)

# # align sample order
# id <- id %>% filter(id$IID %in% cwp$PublicID)
# if(all.equal(cwp$PublicID, id$IID) != TRUE){
#   cwp <- cwp[match(id$IID, cwp$PublicID),]
# }
pheno_val <- cwp[match(id_val, cwp$PublicID),]
all.equal(pheno_val$PublicID, id_val)

covar_val <- covar[match(id_val, covar$PublicID),]
all.equal(covar_val$PublicID, id_val)
colnames(covar_val)[1] <- 'IID'


out_all$pgs$`1`[,1]

v <- validate(out_all, keep=keep_val, covar=covar_val, pheno=pheno_val$case)

id_val_test <- id_snp[!keep_ref]
row_test <- which(!id_val_test %in% c(id_val))
id_test <- id_val_test[row_test]


pheno_test <- cwp[match(id_val_test, cwp$PublicID),]

glm(pheno_test$case~out_all$pgs$`1`[,1], family = "binomial") %>% summary



# the rest are split into two, one as validation, one as test, then reverse, similar to 2-fold cross-validation
# the PRS in test split are calculated and stacked
# 
# # apply to independent test datac(rep(F,1500),rep(T,N-1500))
# out2 <- subset(out, s=v$best.s, lambda=v$best.lambda)
# v2 <- validate(out2, test.bfile=test.bfile, keep = c(rep(F,1500),rep(T,N-1500)), 
#                pheno=cwp$case[-c(1:1500)])
# 
# 
# v <- validate(out)
# 
# glm(cwp$case[-c(1:1500)]~v2$pgs[[1]], family = "binomial") %>% summary