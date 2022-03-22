library(dplyr)
library(magrittr)
library(DescTools)
library(ggplot2)
library(pROC)


cwp <- read.table('id_AB_test_cwp.txt', header = T)
covar <- read.table('Covar.txt', header = T)

# PRS for cwp
prs <- read.table('./CWP_test/outp.all_score', header = T)
prs %<>% select(IID, Pt_5e.08)
colnames(prs) <- c('IID', 'CWP')

# PRS for glycans
files <- list.files(path = './gly_test/', pattern = "IGP.+.all_score")
prs_gly <- lapply(files, function(e) read.table(paste('./gly_test/',e, sep = ''), header = T, 
                                                col.names=c('FID','IID',strsplit(e,'\\.')[[1]][1])))

prs_gly <- do.call(cbind, prs_gly)
IID <- prs_gly$IID
prs_gly %<>% select(starts_with('IGP'))
prs_gly %<>% mutate(IID=IID, .before='IGP49')
prs_gly %>% head

# (S)(P)O2PLS joint components
Tt <- readRDS('Tp_sel_test.rds')
# Ttt <- readRDS('T_sel_test.rds') #O2PLS on selected SNPs
# Tttt <- readRDS('T_test.rds') # SO2PLS on all SNPs
Tt <- as_tibble(Tt)
colnames(Tt) <- c('T1','T2','T3')
T_id <- read.table('psam_test', header = T)
Tt <- bind_cols(T_id, Tt)

cwp <- Reduce(merge, list(cwp, prs, Tt, prs_gly, covar))

# AUC plot function
auc_plot_test <- function(fit_obj){
  prob=predict(fit_obj,type=c("response"))
  g <- roc(fit_obj$y ~ prob)
  ggroc(g, alpha = 1, legacy.axes = TRUE) + xlab("FPR") + ylab("TPR") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed") +
    theme(legend.position = "none") + 
    geom_text(x=0.7, y=0.3, size = 6, label=paste("AUC =", round(g$auc, digit = 3)))
}

# only covariates
fit_cov <- glm(case~Sex+Age+Alc+Smk+Bmi, data = cwp, family = 'binomial')
fit_cov %>% summary
PseudoR2(fit_cov,'Nagelkerke')
auc_plot_test(fit_cov)

# covariates + PRS-cwp
fit_prs <- glm(case~Sex+Age+Alc+Smk+Bmi+CWP, data = cwp, family = 'binomial')
fit_prs %>% summary
PseudoR2(fit_prs,'Nagelkerke')
auc_plot_test(fit_prs)


# covariates + PRS_cwp + PRS_gly
fit_prss <- glm(case~T1+T2+T3, data = cwp, family = 'binomial')
fit_prss %>% summary
PseudoR2(fit_prss,'Nagelkerke')
auc_plot_test(fit_prss)


# covariates + PRS_cwp + jointPC
fit_prso2 <- glm(case~Bmi, data = cwp, family = 'binomial')
fit_prso2 %>% summary
PseudoR2(fit_prso2,'Nagelkerke')
auc_plot_test(fit_prso2)

