library(readxl)
alcohol1 <- readxl::read_excel('/home/z/Data/TwinsUK_CWP/E1124_01062021/E1124_Alcohol.xlsx', 
                               sheet = 1) 
alcohol1$DOV <- as.Date(alcohol1$DOV,'%d/%m/%Y')
alcohol1$DOV <- as.numeric(format(alcohol1$DOV,'%Y'))
alcohol1 <- alcohol1[-1,-c(4,6)]
colnames(alcohol1)[c(3,4)] <- c("current_alc_habit", "life_alc_habit")
table(alcohol1$current_alc_habit, useNA = "ifany")
alcohol1$current_alc_habit <- ifelse(alcohol1$current_alc_habit==0,0,
                                     ifelse(alcohol1$current_alc_habit>=1 & alcohol1$current_alc_habit<=2,2,NA))
table(alcohol1$current_alc_habit, useNA = "ifany") #0=814; 2=7539(current user)
table(alcohol1$life_alc_habit)
alcohol1$life_alc_habit <- ifelse(alcohol1$life_alc_habit==0,0,
                                  ifelse(alcohol1$life_alc_habit>=1 & alcohol1$life_alc_habit<=2,1,NA))
table(alcohol1$life_alc_habit) #0=272; 1=5162

alcohol1$current_alc_use <- rowSums(alcohol1[, c(3,4) ], na.rm=TRUE) * ifelse(
  rowSums(is.na(alcohol1[, c(3,4) ])) == 
    ncol(alcohol1[, c(3,4) ]), NA, 1)
alcohol1$current_alc_use <- ifelse(alcohol1$current_alc_use>=2,2,
                                   ifelse(alcohol1$current_alc_use==1,1,0))
table(alcohol1$current_alc_use) #0=546, 1=274 (ex-user), 2=7539 (current user) 
alcohol1 <- alcohol1[,c(1,2,5)]
alcohol1$current_alc_use <- as.factor(as.character(alcohol1$current_alc_use))
levels(alcohol1$current_alc_use) <- c("never","ex_user", "current_user")



alcohol2 <- readxl::read_excel('/home/z/Data/TwinsUK_CWP/E1124_01062021/E1124_Alcohol.xlsx', sheet = 2)
colnames(alcohol2)[2] <- "DOV" 
alcohol2 <- alcohol2[-1,c(1:4)]

table(alcohol2$Q2_9) #Are you non-drinker/teetotaler now? (i.e. You have not drunk any alcohol for the past three months), 1=yes, 0=no
table(alcohol2$Q2_10) # Have yoy ever drunk alcohol in the past?, 1= yes, 0=no

alcohol2$drinker_now <- ifelse(alcohol2$Q2_9==1,0,2)
table(alcohol2$drinker_now) # 0 = 909 (non-drinker), 2=4981 means drinker now


alcohol2$Q2_10 <- as.numeric(alcohol2$Q2_10)
alcohol2$current_alc_use <- rowSums(alcohol2[,c(4,5) ], na.rm=TRUE) * ifelse(
  rowSums(is.na(alcohol2[, c(4,5) ])) == 
    ncol(alcohol2[, c(4,5) ]), NA, 1)

table(alcohol2$current_alc_use)

alcohol2$current_alc_use <- ifelse(alcohol2$current_alc_use>=2,2,
                                   ifelse(alcohol2$current_alc_use==1,1,0))


alcohol2$current_alc_use <- as.factor(as.character(alcohol2$current_alc_use))                              
levels(alcohol2$current_alc_use) <- c("never","ex_user", "current_user")
alcohol2 <- alcohol2[,c(1,2,6)]



alcohol6 <- readxl::read_excel('/home/z/Data/TwinsUK_CWP/E1124_01062021/E1124_Alcohol.xlsx', sheet = 6)
colnames(alcohol6)[c(1,2)] <- c("Public_ID","DOV") 
alcohol6 <- alcohol6[-1,]
alcohol6[3:6] <- sapply(alcohol6[3:6],as.numeric)
sapply(alcohol6, class)

alcohol6$current_alc_use <- rowSums(alcohol6[, c(3:6) ], na.rm=TRUE) * ifelse(
  rowSums(is.na(alcohol6[, c(3:6) ])) == 
    ncol(alcohol6[, c(3:6) ]), NA, 1)

alcohol6$current_alc_use <- ifelse(alcohol6$current_alc_use>=1,"current/exuser","never")
alcohol6 <- alcohol6[,c(1,2,7)]
table(alcohol6$current_alc_use, useNA = "ifany")




alcohol7 <- readxl::read_excel('/home/z/Data/TwinsUK_CWP/E1124_01062021/E1124_Alcohol.xlsx', sheet = 7)
colnames(alcohol7)[c(1,2)] <- c("Public_ID","DOV") 

alcohol7$DOV <- as.Date(alcohol7$DOV,'%d/%m/%Y')
alcohol7$DOV <- as.numeric(format(alcohol7$DOV,'%Y'))

alcohol7 <- alcohol7[-1,]
alcohol7[3:8] <- sapply(alcohol7[3:8],as.numeric)
sapply(alcohol7, class)
alcohol7$current_alc_use <- rowSums(alcohol7[, c(3:8) ], na.rm=TRUE) * ifelse(
  rowSums(is.na(alcohol7[, c(3:8) ])) == 
    ncol(alcohol7[, c(3:8) ]), NA, 1)
alcohol7$current_alc_use <- ifelse(alcohol7$current_alc_use>=1,"current/exuser","never")
alcohol7 <- alcohol7[,c(1,2,9)]
table(alcohol7$current_alc_use, useNA = "ifany")



alcohol5 <- readxl::read_excel('/home/z/Data/TwinsUK_CWP/E1124_01062021/E1124_Alcohol.xlsx', sheet = 5)
colnames(alcohol5)[c(2,5)] <- c("DOV","Q11A_138")
alcohol5$DOV <- as.Date(alcohol5$DOV,'%d/%m/%Y')
alcohol5$DOV <- as.numeric(format(alcohol5$DOV,'%Y'))
alcohol5 <- alcohol5[-1,c(1:5)]
table(alcohol5$Q11A_136)
alcohol5$Q11A_137 <- ifelse(alcohol5$Q11A_137>=1,1,0)
table(alcohol5$Q11A_137, useNA = "ifany")
table(alcohol5$Q11A_138, useNA="ifany")
alcohol5$Q11A_138 <- ifelse(alcohol5$Q11A_138==0,0,
                            ifelse(alcohol5$Q11A_138>=1&alcohol5$Q11A_138<=6,1,NA))

alcohol5$Q11A_136 <- as.numeric(alcohol5$Q11A_136)

alcohol5$alcohol_use <- rowSums(alcohol5[, c(4,5) ], na.rm=TRUE) * ifelse(
  rowSums(is.na(alcohol5[, c(4,5) ])) == 
    ncol(alcohol5[, c(4,5) ]), NA, 1)
table(alcohol5$alcohol_use, useNA = "ifany")

alcohol5$alcohol_use <- ifelse(alcohol5$alcohol_use>=1,2,0)
table(alcohol5$alcohol_use)

alcohol5$current_alc_use <- rowSums(alcohol5[, c(3,6) ], na.rm=TRUE) * ifelse(
  rowSums(is.na(alcohol5[, c(3,6) ])) == 
    ncol(alcohol5[, c(3,6) ]), NA, 1)
table(alcohol5$current_alc_use)
alcohol5$current_alc_use <- ifelse(alcohol5$current_alc_use>=2,2,
                                   ifelse(alcohol5$current_alc_use==1,1,0))
alcohol5$current_alc_use <- as.factor(as.character(alcohol5$current_alc_use))                              
levels(alcohol5$current_alc_use) <- c("never","ex_user", "current_user")
# never user: 481, exuser: 682, current user: 4742
alcohol5 <- alcohol5[,c(1,2,7)]


colnames(alcohol6)[1] <- "PublicID"
colnames(alcohol7)[1] <- "PublicID"

alc_merged <- rbind(alcohol1,alcohol2,alcohol5,alcohol6,alcohol7)

alc_merged = alc_merged %>%
  mutate(duplicated = duplicated(alc_merged$PublicID, fromLast = T)) %>% 
  mutate(duplicated = ifelse(duplicated== T, "repl_yes", "repl_no"))

table(alc_merged$duplicated)

# dup the old records
alc <- alc_merged[!duplicated(alc_merged$PublicID, fromLast = TRUE), ]
saveRDS(alc, file='/home/z/Data/TwinsUK_CWP/E1124_01062021/alcohol.rds')
