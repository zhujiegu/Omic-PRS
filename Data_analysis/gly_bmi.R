library(dplyr)
library(readxl)
library(magrittr)
library(glmnet)
library(moments)

bmi <- readRDS('bmi.rds')

gly <- readRDS('../gly_gp_orcades.rds')
colnames(gly) <- c('IID', paste0('IGP',1:23))
gly_bmi <- merge(gly, bmi)

# single point
p_single <- c()
for(gp in colnames(gly)[-1]){
  p_single <- append(p_single, (lm(paste0('bmi_log~',gp), data = gly_bmi) %>% summary)$coef[2,4])
}
names(p_single) <- colnames(gly)[-1]
p_single %>% sort

sel <- names(p_single)[p_single<0.05]
gly_sel <- gly[c("IID",sel)] 

saveRDS(gly_sel, file = 'gly_sel_bmi.rds')

#################################################
# PCA
gly_sel_bmi <- merge(gly_sel, bmi, by='IID')
gly_pca <- prcomp(gly_sel_bmi %>% dplyr::select(starts_with('IGP')))
plot(gly_pca$sdev/sum(gly_pca$sdev))
cumsum(gly_pca$sdev)/sum(gly_pca$sdev)

gly_pc <- data.frame(IID=gly_sel_bmi$IID, gly_pca$x[,1:5], bmi_log=gly_sel_bmi$bmi_log)
lm(bmi_log~.-IID, gly_pc) %>% summary()

id_train <- read.table('../id_train')
id_test <- read.table('../id_test')
lm(bmi_log~.-IID, filter(gly_pc, IID%in%id_train$V1)) %>% summary()
lm(bmi_log~.-IID, filter(gly_pc, IID%in%id_test$V1)) %>% summary()
