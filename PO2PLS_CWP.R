library(dplyr)
library(magrittr)

jointPC <- readRDS('jointPC_490.rds')
cwp <- readRDS(file='cwp.rds')

jointPC %<>% mutate(cwp = cwp$case)

all.equal(cwp$PublicID, rownames(jointPC))

glm(cwp~., family = 'binomial', data=jointPC) %>% summary


cor(jointPC)
cor(jointPC, fit_o2$Tt)

glm(cwp$case~fit_o2$Tt, family = 'binomial') %>% summary
