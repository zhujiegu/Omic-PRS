library(dplyr)
library(readxl)

# genetic 5654
###############################################
id <- read.table('chr22.psam')
###############################################


# glycans 4625
###############################################
gly <- read.csv('~/Data/TwinsUK/Glycans_PID/glycans.igg.global.combat.scale.processed_PID.csv', stringsAsFactors = F)
gly_meta <- read_excel('~/Data/TwinsUK/E995_310519_PID.xlsx') %>%
  filter(PublicID %in% gly$PublicID) %>% arrange(PublicID)

# There are duplicates, delete the old one from batch 1(information from Genos)
gly <- gly[-which(gly$PublicID %>% duplicated(fromLast=T)),]
gly %<>% filter(PublicID %in% gly_meta$PublicID) %>% arrange(PublicID)
all.equal(gly$PublicID, gly_meta$PublicID)
###############################################


# CWP 9257
###############################################
cwp <- readRDS('cwp_all.rds')
###############################################


# overlap
###############################################
id %>% filter(V1 %in% cwp$PublicID) %>% nrow # 5447 genetics have CWP
id %>% filter(V1 %in% gly$PublicID) %>% nrow # 4478 genetics have gly
id %>% filter(V1 %in% gly$PublicID) %>% filter(V1 %in% cwp$PublicID) %>% nrow 
# all overlaping genetics and gly have CWP
gly %>% filter(PublicID %in% cwp$PublicID) %>% nrow # all gly have CWP
###############################################

# id lists
# A: SNP, gly, cwp
# B: SNP, cwp
# C: SNP
###############################################
id_A <- (id %>% filter(V1 %in% gly$PublicID))$V1
id_B <- (id %>% filter(!V1 %in% gly$PublicID) %>% filter(V1 %in% cwp$PublicID))$V1
id_C <- (id %>% filter(!V1 %in% gly$PublicID) %>% filter(!V1 %in% cwp$PublicID))$V1
###############################################

# decouple twins
###############################################
sample_famC <- substr(id_C, 1, 10)
sample_famB <- substr(id_B, 1, 10)
sample_famA <- substr(id_A, 1, 10)

# remove twin starting from C, then B, last in A
#C
id_C <- id_C[!sample_famC %in% c(sample_famA, sample_famB)]
sample_famC <- substr(id_C, 1, 10)
sel_twinC <- !duplicated(sample_famC)

#B
id_B <- id_B[!sample_famB %in% sample_famA]
sample_famB <- substr(id_B, 1, 10)
sel_twinB <- !duplicated(sample_famB)

#A
sample_famA <- substr(id_A, 1, 10)
sel_twinA <- !duplicated(sample_famA)

id_A_decpl <- data.frame(IID=id_A[sel_twinA])  #2832
id_B_decpl <- data.frame(IID=id_B[sel_twinB])  #633
id_C_decpl <- data.frame(IID=id_C[sel_twinC])  #136

id_decpl <- bind_rows(id_A_decpl, id_B_decpl, id_C_decpl)
head(id_decpl)
write.table(id_decpl, file='id_decpl.txt', col.names = T, row.names = F, quote = F)
write.table(id_A_decpl, file='id_A_decpl.txt', col.names = T, row.names = F, quote = F)
write.table(id_B_decpl, file='id_B_decpl.txt', col.names = T, row.names = F, quote = F)
write.table(id_C_decpl, file='id_C_decpl.txt', col.names = T, row.names = F, quote = F)

###############################################

cwp %>% filter(PublicID %in% id_overlap_decpl$IID) %>% select(case) %>% sum
