library(glmnet)
library(moments)
library(stringr)
library(dplyr)

# outcome bmi
bmi <- readRDS('bmi.rds')
pheno <- read.table('/home/z/Data/orcades_data_for_zhujie/phenotypes/orcades_phenotypes_jay.tsv', sep = '\t', header = T)

# screened data or not?
screened <- T
################################################################
# cov BMI
lm_cov <- lm(bmi_log~age+sex, data=bmi)
summary(lm_cov)
################################################################
# PRS BMI
bmi_prs <- read.table('./fit_files/PRS_bmi.best', header = T)
bmi_prs <- merge(bmi_prs, bmi, by='IID')

# # best p-value threshold
# p_thresh <- function(dat){
#   thresh <- colnames(dat %>% dplyr::select(starts_with('Pt')))
#   r2_tmp <- c()
#   for (j in thresh) {
#     sumy <- lm(paste0('bmi_log~',j,'+age+sex'), data = dat) %>% summary
#     r2_tmp <- append(r2_tmp, sumy$r.squared)
#   }
#   print(r2_tmp)
#   return(thresh[which.max(r2_tmp)])
# }
# p_best <- p_thresh(bmi_prs)

lm_bmi <- lm(bmi_log~PRS, data=bmi_prs)
# lm_bmi <- lm(bmi_log~Pt_1+age+sex, data=bmi_prs)
summary(lm_bmi)
############################################################################################
# glycomics
############################################################################################
# PRS 
# Construct PRS matrix
gly_sel <- readRDS('gly_sel_bmi.rds')
igp_id <- colnames(gly_sel)[-1]
bmi_gly_l <- vector(mode = 'list', length = length(igp_id))

for(i in 1:length(igp_id)){
  bmi_gly_l[[i]] <- read.table(paste0('./fit_files_gly/PRS_',igp_id[i],'.best'), header = T)
}

gly_prs <- sapply(1:length(igp_id), function(j){
  bmi_gly_l[[j]]$PRS
})
gly_prs <- data.frame(IID=bmi_gly_l[[1]]$IID, IGP=gly_prs)
colnames(gly_prs) <- c('IID',igp_id)
gly_prs %>% dplyr::select(starts_with('IGP')) %>% boxplot
gly_prs %>% dplyr::select(starts_with('IGP')) %>% var %>% diag

# recover variance of each glycan
gly_prs[,-1] <- gly_prs[,-1] %>% scale
# gly_prs[,-1] <- as.matrix(gly_prs[,-1]) %*% diag(apply(gly_sel[,-1],2,sd))
gly_prs %>% dplyr::select(starts_with('IGP')) %>% boxplot
gly_prs <- merge(gly_prs, bmi, by='IID')

lm(bmi_log~.-IID-age-sex, data=gly_prs) %>% summary

# correlation with gly data
gly_sel <- merge(gly_sel, gly_prs, by='IID')
corr_gly <- sapply(igp_id, function(i) cor(gly_sel[c(paste0(i,'.x'), paste0(i,'.y'))])[1,2])
corr_gly %>% barplot()
# plot 
cairo_pdf(file = "/home/z/Dropbox/articles/Omics_PRS/PRS_cor_gly.pdf", width = 8, height =4, onefile =F, fallback_resolution = 600)
corr_gly %>% barplot(ylab='correlation with data')
dev.off()

# PCA
gly_pca <- prcomp(gly_prs %>% dplyr::select(starts_with('IGP')))
plot(gly_pca$sdev/sum(gly_pca$sdev))
cumsum(gly_pca$sdev)/sum(gly_pca$sdev)

gly_pc <- data.frame(gly_pca$x[,1:5], bmi_log=gly_prs$bmi_log)
# gly_pc <- data.frame(gly_pca$x[,1:6], bmi_log=gly_prs$bmi_log, age=bmi$age, sex=bmi$sex)
lm_prs_gly <- lm(bmi_log~., gly_pc) 
lm_prs_gly %>% summary()

################################################################
# O2PLS
if(screened){
  fit_o2_gly <- readRDS('./fit_files_gly/fit_o2_gly_screen.rds')
  fit_o2snp_gly <- readRDS('./fit_files_gly/fit_o2snp_gly_screen.rds')
}else{
  fit_o2_gly <- readRDS('./fit_files_gly/fit_o2_gly.rds')
  fit_o2snp_gly <- readRDS('./fit_files_gly/fit_o2snp_gly.rds')
}

T_o2_gly <- data.frame(IID=rownames(fit_o2_gly$Tt), Tt=fit_o2_gly$Tt)
T_o2_gly <- merge(T_o2_gly, bmi, by='IID')
# lm_o2 <-lm(bmi_log~.-IID, data = T_o2) 
lm_o2_gly <-lm(bmi_log~.-IID-age-sex, data = T_o2_gly) 
lm_o2_gly %>% summary
# lm(bmi_log~.-IID, data = T_o2) %>% summary

# O2PLS snp
T_o2snp_gly <- data.frame(IID=rownames(fit_o2snp_gly$Tt), Tt=fit_o2snp_gly$Tt)
T_o2snp_gly <- merge(T_o2snp_gly, bmi, by='IID')
# lm_o2snp <- lm(bmi_log~.-IID, data = T_o2snp) 
lm_o2snp_gly <-lm(bmi_log~.-IID-age-sex, data = T_o2snp_gly) 
lm_o2snp_gly %>% summary

################################################################
# S/GO2PLS
if(screened){
  fit_so2_gly <- readRDS('./fit_files_gly/fit_so2_gly_screen.rds')
  fit_so2snp_gly <- readRDS('./fit_files_gly/fit_so2snp_gly_screen.rds')
  fit_go2_gly <- readRDS('./fit_files_gly/fit_go2_gly_screen.rds')
}else{
  fit_so2_gly <- readRDS('./fit_files_gly/fit_so2_gly.rds')
  fit_so2snp_gly <- readRDS('./fit_files_gly/fit_so2snp_gly.rds')
  fit_go2_gly <- readRDS('./fit_files_gly/fit_go2_gly.rds')
}

# SO2PLS
T_so2_gly <- data.frame(IID=rownames(fit_so2_gly$Tt), Tt=fit_so2_gly$Tt)
T_so2_gly <- merge(T_so2_gly, bmi, by='IID')
lm_so2_gly <-lm(bmi_log~.-IID-age-sex, data = T_so2_gly) 
lm_so2_gly %>% summary
# lm(bmi_log~.-IID, data = T_so2_gly) %>% summary

# so2PLS snp
T_so2snp_gly <- data.frame(IID=rownames(fit_so2snp_gly$Tt), Tt=fit_so2snp_gly$Tt)
T_so2snp_gly <- merge(T_so2snp_gly, bmi, by='IID')
lm_so2snp_gly <-lm(bmi_log~.-IID-age-sex, data = T_so2snp_gly) 
lm_so2snp_gly %>% summary

# GO2PLS
T_go2_gly <- data.frame(IID=rownames(fit_go2_gly$Tt), Tt=fit_go2_gly$Tt)
T_go2_gly <- merge(T_go2_gly, bmi, by='IID')
lm_go2_gly <-lm(bmi_log~.-IID-age-sex, data = T_go2_gly) 
lm_go2_gly %>% summary
# lm(bmi_log~.-IID, data = T_go2_gly) %>% summary

################################################################
# PO2PLS
if(screened){
  T_p <- readRDS('./fit_files_gly/T_p_6_screen.rds')
}else{
  T_p <- readRDS('./fit_files_gly/T_p_2.rds')
}

T_p <- data.frame(IID=rownames(T_p), Tt=T_p)
T_p <- merge(T_p, bmi, by='IID')
# lm_p <- lm(bmi_log~.-IID, data = T_p)
lm_p_gly <- lm(bmi_log~.-IID-sex-age, data = T_p)
lm_p_gly %>% summary

# lm(bmi_log~.-IID-age-sex, data = T_p) %>% summary


############################################################################################
# Mtb
############################################################################################
mtb <- readRDS('mtb_sel_bmi.rds')

# PRS 
# Construct PRS matrix
fls <- list.files(path = './fit_files_mtb/')
fls %<>% subset(endsWith(fls, ".all_score"))

mtbs <- sapply(fls, function(e) str_sub(e, 5, -11))
bmi_mtb_l <- lapply(colnames(mtb)[-1], 
       function(e) read.table(paste0('./fit_files_mtb/PRS_',e,'.best'), header = T))

mtb_prs <- sapply(bmi_mtb_l, function(e) e$PRS)

mtb_prs <- data.frame(IID=bmi_mtb_l[[1]]$IID, mtb=mtb_prs)
colnames(mtb_prs) <- colnames(mtb)
mtb_prs %>% dplyr::select(-IID) %>% boxplot
# mtb_prs %>% dplyr::select(starts_with('IGP')) %>% var %>% diag

mtb_prs[,-1] <- mtb_prs[,-1] %>% scale
# mtb_prs[,-1] <- as.matrix(mtb_prs[,-1]) %*% diag(apply(mtb_sel[,-1],2,sd))
mtb_prs %>% dplyr::select(-IID) %>% boxplot
mtb_prs <- merge(mtb_prs, bmi, by='IID')

# correlation with metabolomics data
mtb_test <- merge(mtb, mtb_prs, by='IID')
corr_mtb <- sapply(colnames(mtb)[-1], function(i) cor(mtb_test[c(paste0(i,'.x'), paste0(i,'.y'))])[1,2])
corr_mtb %>% barplot(ylab='correlation with data', las=2, cex.axis=2,cex.lab=2)
# plot 
cairo_pdf(file = "/home/z/Dropbox/articles/Omics_PRS/PRS_cor_mtb.pdf", width = 16, height =8, onefile =F, fallback_resolution = 600)
corr_mtb %>% barplot(ylab='correlation with data', las=2, cex.axis=1.5,cex.lab=1.5)
dev.off()
# lm(bmi_log~.-IID-age-sex, mtb_prs) %>% summary

# PCA
mtb_pca <- prcomp(mtb_prs %>% dplyr::select(-c('IID','age','sex', 'bmi_log')))
plot(mtb_pca$sdev/sum(mtb_pca$sdev))
cumsum(mtb_pca$sdev)/sum(mtb_pca$sdev)

mtb_pc <- data.frame(mtb_pca$x[,1:17], bmi_log=mtb_prs$bmi_log)
# mtb_pc <- data.frame(mtb_pca$x[,1:6], bmi_log=mtb_prs$bmi_log, age=bmi$age, sex=bmi$sex)
lm_prs_mtb <- lm(bmi_log~., mtb_pc) 
lm_prs_mtb %>% summary()

################################################################
# O2PLS
if(screened){
  fit_o2_mtb <- readRDS('./fit_files_mtb/fit_o2_mtb_screen.rds')
  fit_o2snp_mtb <- readRDS('./fit_files_mtb/fit_o2snp_mtb_screen.rds')
}else{
  fit_o2_mtb <- readRDS('./fit_files_mtb/fit_o2_mtb.rds')
  fit_o2snp_mtb <- readRDS('./fit_files_mtb/fit_o2snp_mtb.rds')
}

T_o2_mtb <- data.frame(IID=rownames(fit_o2_mtb$Tt), Tt=fit_o2_mtb$Tt)
T_o2_mtb <- merge(T_o2_mtb, bmi, by='IID')
# lm_o2 <-lm(bmi_log~.-IID, data = T_o2) 
lm_o2_mtb <-lm(bmi_log~.-IID-age-sex, data = T_o2_mtb) 
lm_o2_mtb %>% summary

# O2PLS snp
T_o2snp_mtb <- data.frame(IID=rownames(fit_o2snp_mtb$Tt), Tt=fit_o2snp_mtb$Tt)
T_o2snp_mtb <- merge(T_o2snp_mtb, bmi, by='IID')
# lm_o2snp <-lm(bmi_log~.-IID, data = T_o2snp) 
lm_o2snp_mtb <-lm(bmi_log~.-IID-age-sex, data = T_o2snp_mtb) 
lm_o2snp_mtb %>% summary

################################################################
# S/GO2PLS
if(screened){
  fit_so2_mtb <- readRDS('./fit_files_mtb/fit_so2_mtb_screen.rds')
  fit_so2snp_mtb <- readRDS('./fit_files_mtb/fit_so2snp_mtb_screen.rds')
  fit_go2_mtb <- readRDS('./fit_files_mtb/fit_go2_mtb_screen.rds')
}else{
  fit_so2_mtb <- readRDS('./fit_files_mtb/fit_so2_mtb.rds')
  fit_so2snp_mtb <- readRDS('./fit_files_mtb/fit_so2snp_mtb.rds')
  fit_go2_mtb <- readRDS('./fit_files_mtb/fit_go2_mtb.rds')
}

# SO2PLS
T_so2_mtb <- data.frame(IID=rownames(fit_so2_mtb$Tt), Tt=fit_so2_mtb$Tt)
T_so2_mtb <- merge(T_so2_mtb, bmi, by='IID')
lm_so2_mtb <-lm(bmi_log~.-IID-age-sex, data = T_so2_mtb) 
lm_so2_mtb %>% summary
# lm(bmi_log~.-IID, data = T_so2_mtb) %>% summary

# so2PLS snp
T_so2snp_mtb <- data.frame(IID=rownames(fit_so2snp_mtb$Tt), Tt=fit_so2snp_mtb$Tt)
T_so2snp_mtb <- merge(T_so2snp_mtb, bmi, by='IID')
lm_so2snp_mtb <-lm(bmi_log~.-IID-age-sex, data = T_so2snp_mtb) 
lm_so2snp_mtb %>% summary

# GO2PLS
T_go2_mtb <- data.frame(IID=rownames(fit_go2_mtb$Tt), Tt=fit_go2_mtb$Tt)
T_go2_mtb <- merge(T_go2_mtb, bmi, by='IID')
lm_go2_mtb <-lm(bmi_log~.-IID-age-sex, data = T_go2_mtb) 
lm_go2_mtb %>% summary
# lm(bmi_log~.-IID, data = T_go2_mtb) %>% summary


################################################################
# PO2PLS
if(screened){
  T_p_mtb <- readRDS('./fit_files_mtb/T_p_5_screen.rds')
}else{
  T_p_mtb <- readRDS('./fit_files_mtb/T_p_10.rds')
}

T_p_mtb <- data.frame(IID=rownames(T_p_mtb), Tt=T_p_mtb)
T_p_mtb <- merge(T_p_mtb, bmi, by='IID')
# lm_p <- lm(bmi_log~.-IID, data = T_p)
lm_p_mtb <- lm(bmi_log~.-IID-sex-age, data = T_p_mtb)
lm_p_mtb %>% summary

lm(bmi_log~Tt.1, data = T_p_mtb)%>% summary
# lm_p_mtb %>% summary

qq_p_mtb <- T_p_mtb %>% arrange(Tt.1)
t.test(head(qq_p_mtb, 150)$bmi_log, tail(qq_p_mtb, 150)$bmi_log)
################################################################
# # combine
# cmb <- merge(T_o2, bmi_prs, by='IID')
# lm(bmi_log.x~Tt.1+Tt.2+Tt.3+Tt.4+Pt_1, data = cmb) %>% summary
# 
# cmb <- merge(T_p, bmi_prs, by='IID')
# lm(bmi_log.x~Tt.1+Tt.2+Tt.3+Tt.4+Pt_1, data = cmb) %>% summary
# 
# 
# cmb <- cbind(gly_pc, bmi_prs$Pt_1)
# lm(bmi_log~., data = cmb) %>% summary

################################################################
# save regression coefficient for TwinsUK
save(lm_prs_gly, lm_o2_gly, lm_o2snp_gly, lm_so2_gly, lm_so2snp_gly, lm_go2_gly, lm_p_gly, 
     lm_prs_mtb, lm_o2_mtb, lm_o2snp_mtb, lm_so2_mtb, lm_so2snp_mtb, lm_go2_mtb, lm_p_mtb,
     gly_pca, mtb_pca, file = 'lmfit_orcades.RData')

# save(lm_o2_gly, lm_o2snp_gly, lm_p_gly, lm_o2_mtb, lm_o2snp_mtb, lm_p_mtb, file = 'lmfit_orcades_screen.RData')
