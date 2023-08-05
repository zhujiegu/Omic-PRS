library(dplyr)
library(moments)
library(magrittr)
library(ggplot2)

id_norel <- read.table('../chr_all_QC_25.king.cutoff.in.id')

pheno <- read.table('/home/z/Data/orcades_data_for_zhujie/phenotypes/orcades_phenotypes_jay.tsv', sep = '\t', header = T)
pheno %<>% select(iid, bmi,age,sex)
pheno <- pheno %>% filter(iid %in% id_norel$V1)
pheno %<>% filter(!is.na(bmi)) %>% filter(!is.na(age)) %>% filter(!is.na(sex))
# check outlier
pheno$bmi %>% density %>% plot
# skewness
pheno$bmi %>% skewness()
pheno$bmi %>% log %>% skewness()

pheno %<>% mutate(bmi_log = log(bmi)) %>% select(-bmi)
pheno$sex %<>% as.factor()
pheno %>% head
colnames(pheno)[1] <- 'IID'
pheno$bmi_log %>% hist(title())
saveRDS(pheno, file = 'bmi.rds')

cairo_pdf(file = "/home/z/Dropbox/articles/Omics_PRS/BMI_hist_ORCADES.pdf", width = 8, height =4, onefile =F, fallback_resolution = 600)
ggplot(pheno, aes(x=bmi_log)) + 
  geom_histogram(binwidth=0.1, colour="black", fill="white") + 
  theme_bw()+
  ggtitle("Histogram of log(BMI) in ORCADES") +
  scale_x_continuous(breaks = seq(2.9, 4.0, by = 0.1)) + #limits = c(2.6,4.1)
  theme(title = element_text(size=16,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text=element_text(size=14,face="bold"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=12,face="bold"))
dev.off()

# for PRSice-2
# write.table(pheno, file='bmi.txt', col.names = T, row.names = F, quote = F)
