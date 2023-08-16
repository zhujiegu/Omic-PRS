library(dplyr)
library(readxl)
library(magrittr)
library(glmnet)
library(moments)

bmi <- readRDS('bmi.rds')

mtb <- readRDS('../mtb_processed_orcades.rds')
mtb_bmi <- merge(mtb, bmi, by='IID')

# single point
p_single <- c()
for(gp in colnames(mtb)[-1]){
  p_single <- append(p_single, (lm(paste0('bmi_log~',gp), data = mtb_bmi) %>% summary)$coef[2,4])
}
names(p_single) <- colnames(mtb)[-1]
p_single %>% sort

sel <- names(p_single)[p_single<0.05]
mtb_sel <- mtb[c("IID",sel)] 

saveRDS(mtb_sel, file = 'mtb_sel_bmi.rds')
#################################################
# PCA
mtb_sel_bmi <- merge(mtb_sel, bmi, by='IID')
mtb_pca <- prcomp(mtb_sel_bmi %>% dplyr::select(-c('IID','age','sex','bmi_log')))
plot(mtb_pca$sdev/sum(mtb_pca$sdev))
cumsum(mtb_pca$sdev)/sum(mtb_pca$sdev)

# BMI
mtb_pc <- data.frame(IID=mtb_sel_bmi$IID, mtb_pca$x[,1:5], bmi_log=mtb_sel_bmi$bmi_log)
lm(bmi_log~.-IID, mtb_pc) %>% summary()

id_train <- read.table('../id_train')
id_test <- read.table('../id_test')
lm(bmi_log~.-IID, filter(mtb_pc, IID%in%id_train$V1)) %>% summary()
lm(bmi_log~.-IID, filter(mtb_pc, IID%in%id_test$V1)) %>% summary()
