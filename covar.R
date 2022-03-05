library(magrittr)
library(dplyr)
library(readxl)

# samples with genetic data
id <- read.table('chr1.fam')
id <- id$V2

# age and sex (at 2014, the last questionnaire time)
age_sex <- read_excel('/home/z/Data/TwinsUK_CWP/E1124_01062021/E1124_TwinDetails.xlsx', 
                      sheet = 1) %>% filter(PublicID %in% id) %>% arrange(PublicID)
age_sex %<>% mutate(age = 2014-YEAR_BIRTH)

# check duplicated samples
age_sex$PublicID %>% duplicated %>% sum

age_sex$YEAR_BIRTH %>% is.na %>% sum
age_sex$SEX %>% is.na %>% sum

head(age_sex)

Covar <- select(age_sex, PublicID, SEX, age)

# Alcohol
alc <- readRDS('/home/z/Data/TwinsUK_CWP/E1124_01062021/alcohol.rds') %>% 
  filter(PublicID %in% id) %>% arrange(PublicID)

# Smoking
smoke <- readRDS('/home/z/Data/TwinsUK_CWP/E1124_01062021/smoking.rds') %>% 
  filter(PublicID %in% id) %>% arrange(PublicID)

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

weight <- weight %>% filter(PublicID %in% height$PublicID) %>% filter(PublicID %in% id)
height <- height %>% filter(PublicID %in% weight$PublicID) %>% filter(PublicID %in% id)
all.equal(weight$PublicID, height$PublicID)
bmi <- weight$Weight/height$Height^2 *10000
weight <- weight %>% mutate(Bmi = bmi)

## merge
Covar <- merge(Covar, select(alc, PublicID, current_alc_use), by='PublicID', all.x=T)
Covar <- merge(Covar, select(smoke, PublicID, sumstat), by='PublicID', all.x=T)
Covar <- merge(Covar, select(weight, PublicID, Bmi), by='PublicID', all.x=T)

colnames(Covar) <- c('IID','Sex','Age','Alc','Smk','Bmi')
head(Covar)

is.na(Covar$Alc) %>% sum
is.na(Covar$Smk) %>% sum
is.na(Covar$Bmi) %>% sum

levels(Covar$Alc) <- c(levels(Covar$Alc), 'unknown')
levels(Covar$Smk) <- c(levels(Covar$Smk), 'unknown')

Covar$Alc[which(is.na(Covar$Alc))] <- 'unknown'
Covar$Smk[which(is.na(Covar$Smk))] <- 'unknown'
Covar$Bmi[which(is.na(Covar$Bmi))] <- median(Covar$Bmi, na.rm = T)

saveRDS(Covar, file = 'Covar.Rds')
write.table(Covar, file='Covar.txt', col.names = T, row.names = F, quote = F)

id_A_training_cwp$IID %in% Covar$IID 
