library(dplyr)
library(magrittr)
library(MASS)
library(moments)
library(mice)
library(readxl)

mtb <- read.table(file = '/home/z/Data/orcades_data_for_zhujie/nmr_metabolomics/d001_orcades_nmr_phenotypes.tsv.gz', sep = '\t', header = TRUE)
# remove column with no variance
mtb <- mtb[,-2]

# check skewness
apply(mtb[,-1], 2, function(e) skewness(e, na.rm = T)) %>% summary

####################################################################################
# Hagenbeek2020
# missing rate > 5%
mis_pct <- apply(is.na(mtb[,-1]), 2, mean)
mtb %<>% dplyr::select(-c(names(mis_pct[mis_pct>0.05])))
sum(is.na(mtb))
######################
# Zeros
sum(!is.na(mtb)&mtb==0)
zero_pct <- apply((mtb[,-1]==0), 2, mean)
# zero_pct %>% sort(decreasing = T)
# replace zeros with 1/2 of the lowest observation
for(j in 2:ncol(mtb)){
  if(sum(mtb[,j]==0, na.rm = T)>0){
    mtb[!is.na(mtb[,j])&mtb[,j]==0,j] <- min(mtb[!mtb[,j]%in%c(0, NA),j])
  }
}
sum(!is.na(mtb)&mtb==0)
######################
# Boxcox
mtb[,-1] %>% boxplot
apply(mtb[,-1], 2, function(e) skewness(e, na.rm = T)) %>% summary
mtb_25 <- mtb[,-1]^0.25
mtb_25 %>% boxplot
apply(mtb_25, 2, function(e) skewness(e, na.rm = T)) %>% density %>% plot
apply(mtb_25, 2, function(e) var(e, na.rm = T)) %>% density %>% plot

######################
# outlier 6 sigma (Macdonald-Dunlop2022)
sum(is.na(mtb_25))
for(j in 1:ncol(mtb_25)){
  m <- mean(mtb_25[,j], na.rm=T)
  sig <- sd(mtb_25[,j], na.rm=T)
  mtb_25[which(mtb_25[,j]<m-6*sig | mtb_25[,j]>m+6*sig) ,j] <- NA
}
sum(is.na(mtb_25))
######################
# missing
md.pattern(mtb_25)
imp <- mice(mtb_25,m=1,meth='norm.predict',seed=500)
mtb_imp <- complete(imp)
sum(is.na(mtb_imp))
mtb_imp <- bind_cols(iid=mtb$iid, mtb_imp)
######################
# sex and age
pheno <- read.table('/home/z/Data/orcades_data_for_zhujie/phenotypes/orcades_phenotypes_jay.tsv', sep = '\t', header = T)
covar <- pheno %>% dplyr::select(iid, age, sex)
covar$sex %<>% as.factor
mtb_imp <- merge(mtb_imp, covar, by='iid')
mtb_imp$age %>% is.na %>% sum
mtb_imp$sex %>% is.na %>% sum

mtb_imp$age %>% density %>% plot
mtb_imp$sex %>% table

correct <- lapply(colnames(mtb_imp)[2:211], function(j){
  fit <- lm(paste0(j,'~age+sex'), data = mtb_imp)
  return(residuals(fit))
})

correct_matrix <- do.call(cbind, correct)
mtb_imp[,2:211] <- correct_matrix
mtb_imp <- mtb_imp %>% dplyr::select(-age,-sex)
####################################################################################
# overlap with Kettunen2016
mtb_orc <- read.table('/home/z/Data/orcades_data_for_zhujie/nmr_metabolomics/d000_nmr_description.csv', header = T, sep = ',')
mtb_Ket <- read_excel('/home/z/Data/orcades_data_for_zhujie/nmr_metabolomics/Kettunen2016_mtb_list.xlsx')

mtb_Ket$variable <- mtb_Ket$Abbreviation %>% tolower()
mtb_Ket$variable <- gsub('\\.', '_', mtb_Ket$variable)

mtb_com <- merge(mtb_orc, mtb_Ket, by='variable')
# # manual check rest of the variables in Kettunen
# check_list <- mtb_Ket %>% filter(!Abbreviation %in% mtb_com$Abbreviation)

mtb_imp %<>% dplyr::select(iid, any_of(mtb_com$variable))
# change metabolites names to Kettunen ones
colnames(mtb_imp) <- c('IID', mtb_com$Abbreviation[match(colnames(mtb_imp)[-1], mtb_com$variable)])

# scale
mtb_imp[,-1] <- scale(mtb_imp[,-1])
mtb_imp[,-1] %>% boxplot
saveRDS(mtb_imp, file='mtb_processed_orcades.rds')
