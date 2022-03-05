library(dplyr)
library(fmsb)

# read PRS in training
score_train <- read.table('./PRSice_training_no_regres/outp.all_score', header = T)
head(score_train)

# remove FID column
score_train <- score_train[,-1]

# check correlation
ggcorrplot::ggcorrplot(cor(score_train[,-1]))

# read cwp
cwp <- readRDS('cwp_all.rds')
colnames(cwp)[1] <- "IID"

# read covariates
Covar <- readRDS('Covar.Rds')

score_train <- merge(score_train, Covar, by = "IID")
score_train <- merge(score_train, cwp, by = "IID")

head(score_train)

# Null model 
null.model <- glm(case~ Sex+Age+Alc+Smk+Bmi, data = score_train, family = 'binomial') %>% summary
null.r2 <- NagelkerkeR2(null.model)

test.model <- glm(case~ Pt_1 + Sex+Age+Alc+Smk+Bmi, data = score_train) %>%
  summary
(result <- data.table::data.table(
  auto = test.model$r.squared - null.r2,
  null = null.r2
))



patients <- score_train %>% filter(case == 0) %>% select(starts_with('Pt'))
(patients >0) %>% sum
# try summary stats from other paper
