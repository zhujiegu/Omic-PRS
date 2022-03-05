library(magrittr)
library(dplyr)
library(tidyverse)
library(ggcorrplot)
library(readxl)


gly <- read.csv('~/Data/TwinsUK/Glycans_PID/glycans.igg.global.combat.scale.processed_PID.csv', stringsAsFactors = F)
gly_meta <- read_excel('~/Data/TwinsUK/E995_310519_PID.xlsx') %>%
  filter(PublicID %in% gly$PublicID) %>% arrange(PublicID)

# There are duplicates, delete the old one from batch 1(information from Genos)
gly <- gly[-which(gly$PublicID %>% duplicated(fromLast=T)),]
gly %<>% filter(PublicID %in% gly_meta$PublicID) %>% arrange(PublicID)
all.equal(gly$PublicID, gly_meta$PublicID)

# age calculate
date <- gly$date %>% str_split('/') 
dat_num <- sapply(date, function(e) as.numeric(e[[3]])+as.numeric(e[[2]])/12+as.numeric(e[[1]])/365)
age <- round(dat_num - gly_meta$YEAR_BIRTH - 0.5, digits = 1)
gly %<>% mutate(Age = age, Sex = as.factor(gly_meta$SEX))


# Alcohol
alc <- readRDS('/home/z/Data/TwinsUK_CWP/E1124_01062021/alcohol.rds') %>% 
  filter(PublicID %in% gly$PublicID) %>% arrange(PublicID)

# Smoking
smoke <- readRDS('/home/z/Data/TwinsUK_CWP/E1124_01062021/smoking.rds') %>% 
  filter(PublicID %in% gly$PublicID) %>% arrange(PublicID)

# BMI
height <- read_excel("/home/z/Data/TwinsUK_CWP/E1124_01062021/E1124_HeightWeight.xlsx", 
                      sheet = 1)
weight <- read_excel("/home/z/Data/TwinsUK_CWP/E1124_01062021/E1124_HeightWeight.xlsx", 
                     sheet = 2)

height = height %>%
  mutate(duplicated = duplicated(height$PublicID, fromLast = T)) %>% 
  mutate(duplicated = ifelse(duplicated== T, "repl_yes", "repl_no"))
height <- height[!duplicated(height$PublicID, fromLast = TRUE), ]

weight = weight %>%
  mutate(duplicated = duplicated(weight$PublicID, fromLast = T)) %>% 
  mutate(duplicated = ifelse(duplicated== T, "repl_yes", "repl_no"))
weight <- weight[!duplicated(weight$PublicID, fromLast = TRUE), ]

weight <- weight %>% filter(PublicID %in% height$PublicID) %>% filter(PublicID %in% gly$PublicID)
height <- height %>% filter(PublicID %in% weight$PublicID) %>% filter(PublicID %in% gly$PublicID)
all.equal(weight$PublicID, height$PublicID)
bmi <- weight$Weight/height$Height^2 *10000
weight <- weight %>% mutate(Bmi = bmi)

## merge
gly <- merge(gly, select(alc, PublicID, current_alc_use), by='PublicID', all.x=T)
gly <- merge(gly, select(smoke, PublicID, sumstat), by='PublicID', all.x=T)
gly <- merge(gly, select(weight, PublicID, Bmi), by='PublicID', all.x=T)

colnames(gly)[which(names(gly) == "current_alc_use")] = 'Alc'
colnames(gly)[which(names(gly) == "sumstat")] = 'Smk'

is.na(gly$Alc) %>% sum
is.na(gly$Smk) %>% sum
is.na(gly$Bmi) %>% sum

levels(gly$Alc) <- c(levels(gly$Alc), 'unknown')
levels(gly$Smk) <- c(levels(gly$Smk), 'unknown')

gly$Alc[which(is.na(gly$Alc))] <- 'unknown'
gly$Smk[which(is.na(gly$Smk))] <- 'unknown'
gly$Bmi[which(is.na(gly$Bmi))] <- median(gly$Bmi, na.rm = T)

glycan <- gly %>% select(contains('GP'))

# Correct for covariates
correct <- lapply(1:dim(glycan)[2], function(i){
  y <- glycan[,i]
  fit <- lm(y~Age + Sex + Alc + Smk + Bmi, data = gly, na.action = na.exclude)
  return(residuals(fit))
})
glycan <- do.call(cbind, correct)

# Asign names
row.names(glycan) <- gly$PublicID
colnames(glycan) <- gly %>% select(contains('GP')) %>% colnames()

ggcorrplot(cor(glycan))
boxplot(glycan)

covar <- gly %>% select(PublicID,Age,Sex,Alc,Smk,Bmi)

saveRDS(covar, file = 'Covar.Rds')
saveRDS(glycan, file = '/home/z/Data/TwinsUK_CWP/Glycan.Rds')
######################################################
gly <- readRDS('/home/z/Data/TwinsUK_CWP/Glycan.Rds')
gly <- data.frame(gly)

gly_IGP49 <- data.frame(IID=row.names(gly), IGP49=gly$IGP49)
head(gly_IGP49)
write.table(gly_IGP49, '/home/z/Data/IgG_glycans_GWAS/IGP49.txt', col.names = T, row.names = F, quote = F)
