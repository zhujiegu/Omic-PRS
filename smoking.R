setwd("/home/z/Data/TwinsUK_CWP/E1124_01062021")
library("readxl")
#install.packages("lubridate")
library("lubridate")
library(dplyr)

# smoking
smoking1 <- read_excel("E1124_Smoking.xlsx", sheet = 1)
smoking2 <- read_excel("E1124_Smoking.xlsx", sheet = 2)


smoking1 <- smoking1[-1,]
smoking2 <- smoking2[-1,]

colnames(smoking1)[c(2,5)] <- c("DOV", "sumstat")
colnames(smoking2)[c(2,9)] <- c("DOV", "sumstat")


smoking1$DOV <- as.Date(smoking1$DOV,'%d/%m/%Y')
smoking2$DOV <- as.Date(smoking2$DOV,'%d/%m/%Y')

smoking1$DOV <- as.numeric(format(smoking1$DOV,'%Y'))
smoking2$DOV <- as.numeric(format(smoking2$DOV,'%Y'))

smoking1 <- smoking1[,c(1,2,5)]
smoking2 <- smoking2[,c(1,2,9)]


table(smoking2$sumstat, useNA = "ifany")
table(smoking1$sumstat, useNA = "ifany")

smoking1$sumstat <- ifelse(smoking1$sumstat==9,1,
                           ifelse(smoking1$sumstat==8,2,0))

smoking2$sumstat <- ifelse(smoking2$sumstat==7,1,
                           ifelse(smoking2$sumstat==6,2,0))

smoking1$sumstat <- as.factor(as.character(smoking1$sumstat))
smoking2$sumstat <- as.factor(as.character(smoking2$sumstat))

levels(smoking1$sumstat) <- c("non_smoker","ex_smoker","smoker")
levels(smoking2$sumstat) <- c("non_smoker","ex_smoker","smoker")

#smoking1 = 1992 to 1996
#smoking2 = 1996 to 2001






smoking3 <- read_excel("E1124_Smoking.xlsx", sheet = 3)
colnames(smoking3)[c(2,4)] <- c("DOV", "sumstat")
table(smoking3$sumstat, useNA = "ifany")
smoking3 <- smoking3[,c(1,2,4)]
table(smoking3$sumstat)
smoking3$sumstat <- ifelse(smoking3$sumstat==7,1,
                           ifelse(smoking3$sumstat==6,2,0))

smoking3$sumstat <- as.factor(as.character(smoking3$sumstat))
levels(smoking3$sumstat) <- c("non_smoker","ex_smoker","smoker")






smoking4 <- read_excel("E1124_Smoking.xlsx", sheet = 4)
colnames(smoking4)[c(2,7)] <- c("DOV", "sumstat")
table(smoking4$sumstat, useNA = "ifany")
smoking4 <- smoking4[-1,c(1,2,7)]
table(smoking4$sumstat)
smoking4$sumstat <- ifelse(smoking4$sumstat==7,1,
                           ifelse(smoking4$sumstat==6,2,0))

smoking4$sumstat <- as.factor(as.character(smoking4$sumstat))
levels(smoking4$sumstat) <- c("non_smoker","ex_smoker","smoker")






smoking5 <- read_excel("E1124_Smoking.xlsx", sheet = 5)
colnames(smoking5)[c(2,3)] <- c("DOV", "sumstat") #Have you smoked cig at least 1 cig daily in last 30 days?
table(smoking5$sumstat, useNA = "ifany")
smoking5 <- smoking5[-1,c(1,2,3)]
table(smoking5$sumstat)
smoking5$sumstat <- ifelse(smoking5$sumstat==2,1,
                           ifelse(smoking5$sumstat==1,2,0))

smoking5$sumstat <- as.factor(as.character(smoking5$sumstat))
levels(smoking5$sumstat) <- c("non_smoker","ex_smoker","smoker")





smoking6 <- read_excel("E1124_Smoking.xlsx", sheet = 6)
colnames(smoking6)[c(2,3)] <- c("DOV", "sumstat")#Have you smoked cig at least 1 cig daily in last 30 days?
table(smoking6$sumstat, useNA = "ifany")
smoking6 <- smoking6[-1,c(1,2,3)]
table(smoking6$sumstat)
smoking6$sumstat <- ifelse(smoking6$sumstat==2,1,
                           ifelse(smoking6$sumstat==1,2,0))

smoking6$sumstat <- as.factor(as.character(smoking6$sumstat))
levels(smoking6$sumstat) <- c("non_smoker","ex_smoker","smoker")


smoking7 <- read_excel("E1124_Smoking.xlsx", sheet = 7)
smoking7 <- smoking7[-1,c(1,2,8,12)]

# currently smoking daily
smoking7$Q9_71 <- ifelse(smoking7$Q9_71==1,1,
                         ifelse(smoking7$Q9_71==0,0,NA))

# have you ever smoked cigarrets?
smoking7$Q9_67 <- ifelse(smoking7$Q9_67==0,0,
                         ifelse(smoking7$Q9_67==1,1,NA))


table(smoking7$Q9_71,smoking7$Q9_67)

smoking7$sumstat <- ifelse(smoking7$Q9_71==1,2,
                           ifelse(smoking7$Q9_67==1&smoking7$Q9_71!=1,1,
                                  ifelse(smoking7$Q9_67==0&smoking7$Q9_71==0,0,NA)))

table(smoking7$sumstat)
smoking7$sumstat <- as.factor(as.character(smoking7$sumstat))
levels(smoking7$sumstat) <- c("non_smoker","ex_smoker","smoker")

smoking7 <- smoking7[,c(1,2,5)]
colnames(smoking7)[c(2)] <- c("DOV")




smoking8 <- read_excel("E1124_Smoking.xlsx", sheet = 8)
smoking8 <- smoking8[-1,c(1,2,7,9)]
table(smoking8$Q11A_147)
table(smoking8$Q11A_149)

# currently smoking daily
smoking8$Q11A_147 <- ifelse(smoking8$Q11A_147>=1,1,
                            ifelse(smoking8$Q11A_147==0,0,NA))

# have you ever smoked cigarrets?
smoking8$Q11A_149 <- ifelse(smoking8$Q11A_149==0,0,
                            ifelse(smoking8$Q11A_149==1,1,NA))

table(smoking8$Q11A_147,smoking8$Q11A_149)

smoking8$sumstat <- ifelse(smoking8$Q11A_147==1,2,
                           ifelse(smoking8$Q11A_149==1&smoking8$Q11A_147!=1,1,
                                  ifelse(smoking8$Q11A_149==0&smoking8$Q11A_147==0,0,NA)))

table(smoking8$sumstat)
colnames(smoking8)[c(2)] <- c("DOV")
smoking8$DOV <- as.Date(smoking8$DOV,'%d/%m/%Y')
smoking8$DOV <- as.numeric(format(smoking8$DOV,'%Y'))
smoking8 <- smoking8[,c(1,2,5)]
table(smoking8$sumstat)
smoking8$sumstat <- as.factor(as.character(smoking8$sumstat))
levels(smoking8$sumstat) <- c("non_smoker","ex_smoker","smoker")



#2008
smoking9 <- read_excel("E1124_Smoking.xlsx", sheet = 9)
colnames(smoking9)[c(2,3)] <- c("DOV", "sumstat") #Have you smoked cig at least 1 cig daily in last 30 days?
table(smoking9$sumstat, useNA = "ifany")
smoking9 <- smoking9[-1,c(1,2,3)]
smoking9$sumstat <- ifelse(smoking9$sumstat==2,1,
                           ifelse(smoking9$sumstat==1,2,0))
smoking9$sumstat <- as.factor(as.character(smoking9$sumstat))
levels(smoking9$sumstat) <- c("non_smoker","ex_smoker","smoker")


# 2007 to 2010
smoking10 <- read_excel("E1124_Smoking.xlsx", sheet = 10)
smoking10 <- smoking10[-1,c(1,2,3,4)]
# have you ever smoked
table(smoking10$Q17D_11)
smoking10$Q17D_11 <- ifelse(smoking10$Q17D_11==1,1,
                            ifelse(smoking10$Q17D_11==0,0,NA))

# do you currently smoke cigarrates  
table(smoking10$Q17D_12)
smoking10$Q17D_12 <- ifelse(smoking10$Q17D_12==1,1,
                            ifelse(smoking10$Q17D_12==0,0,NA))


smoking10$sumstat <- ifelse(smoking10$Q17D_12==1,2,
                            ifelse(smoking10$Q17D_11==1&smoking10$Q17D_12!=1,1,
                                   ifelse(smoking10$Q17D_11==0&smoking10$Q17D_12==0,0,NA)))


table(smoking10$Q17D_11,smoking10$Q17D_12, useNA = "ifany")

table(smoking10$sumstat)
smoking10$sumstat <- as.factor(as.character(smoking10$sumstat))
levels(smoking10$sumstat) <- c("non_smoker","ex_smoker","smoker")

colnames(smoking10)[c(2)] <- c("DOV")
smoking10$DOV <- as.Date(smoking10$DOV,'%d/%m/%Y')
smoking10$DOV <- as.numeric(format(smoking10$DOV,'%Y'))
smoking10 <- smoking10[,c(1,2,5)]

#2013
smoking11 <- read_excel("E1124_Smoking.xlsx", sheet = 12, skip = 2)

smoking11 <- smoking11[-1,-5]
smoking11$Response_Date <- "2013"

# Have you ever smoked
table(smoking11$Q23_109)
smoking11$Q23_109<- ifelse(smoking11$Q23_109==1,1,
                           ifelse(smoking11$Q23_109==0,0,NA))

# do you currently smoke ciga
smoking11$Q23_110 <- ifelse(smoking11$Q23_110==1,1,
                            ifelse(smoking11$Q23_110==0,0,NA))

table(smoking11$Q23_109,smoking11$Q23_110)
smoking11$sumstat <- ifelse(smoking11$Q23_110==1,2,
                            ifelse(smoking11$Q23_109==1&smoking11$Q23_110!=1,1,
                                   ifelse(smoking11$Q23_109==0&smoking11$Q23_110==0,0,NA)))

table(smoking11$sumstat)
colnames(smoking11)[c(2)] <- c("DOV")
smoking11$sumstat <- as.factor(as.character(smoking11$sumstat))
levels(smoking11$sumstat) <- c("non_smoker","ex_smoker","smoker")
smoking11 <- smoking11[,c(1,2,5)]





colnames(smoking1)[c(1,2,3)] <- c("PublicID","DOV","sumstat")
colnames(smoking2)[c(1,2,3)] <- c("PublicID","DOV","sumstat")
colnames(smoking3)[c(1,2,3)] <- c("PublicID","DOV","sumstat")
colnames(smoking4)[c(1,2,3)] <- c("PublicID","DOV","sumstat")
colnames(smoking5)[c(1,2,3)] <- c("PublicID","DOV","sumstat")
colnames(smoking6)[c(1,2,3)] <- c("PublicID","DOV","sumstat")
colnames(smoking7)[c(1,2,3)] <- c("PublicID","DOV","sumstat")
colnames(smoking8)[c(1,2,3)] <- c("PublicID","DOV","sumstat")
colnames(smoking9)[c(1,2,3)] <- c("PublicID","DOV","sumstat")
colnames(smoking10)[c(1,2,3)] <- c("PublicID","DOV","sumstat")
colnames(smoking11)[c(1,2,3)] <- c("PublicID","DOV","sumstat")



## Merging all smoking data:
smoking_merged <- do.call("rbind", list(smoking1, smoking2, smoking3, smoking4,smoking5,smoking6,
                                        smoking7,smoking8,smoking9,smoking10,smoking11))

table(smoking_merged$sumstat, useNA = "ifany")

smoking_merged = smoking_merged %>%
  mutate(duplicated = duplicated(smoking_merged$PublicID, fromLast = T)) %>% 
  mutate(duplicated = ifelse(duplicated== T, "repl_yes", "repl_no"))

table(smoking_merged$duplicated)

# dup the old records
smoking <- smoking_merged[!duplicated(smoking_merged$PublicID, fromLast = TRUE), ]
saveRDS(smoking, file='/home/z/Data/TwinsUK_CWP/E1124_01062021/smoking.rds')
