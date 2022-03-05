library(readxl)
library(dplyr)

cwp <- read_excel("/home/z/Data/TwinsUK_CWP/E1124_01062021/E1124_pain_q7+11a+14+23+25_lbpdm1.xlsx", sheet =1)
cwp <- cwp[-1,]
cwp <- cwp %>% select(PublicID, case)

saveRDS(cwp, file='cwp_all.rds')
