library(dplyr)
library(magrittr)
library(OmicsPLS)
library(ggplot2)
library(gridExtra)

# load pre-calculated PCA
pca_fit <- readRDS('PCA_snp.rds')
info <- read.table('info_maf.txt', header = T)

cut_pct=0.001

w_l <- vector(mode = "list", length = 3)
for(k in 1:3){
  comp <- (k*5-4):(k*5)
  w <- pca_fit$rotation[,comp]
  for(j in 1:5){
    w[abs(w[,j])<quantile(abs(w[,j]),(1-cut_pct)),j] <- 0
    if(j>1){
      w[,j] <- orth_vec(w[,j], cbind(w[,1:(j-1)]))
    }
    w[,j] <- w[,j]/norm_vec(w[,j])
  }
  w_l[[k]]<-w
}
w <- Reduce(cbind, w_l)


w <- pca_fit$rotation[,1:15]

  for(j in 1:15){
    w[abs(w[,j])<quantile(abs(w[,j]),(1-cut_pct)),j] <- 0
    if(j>1){
      w[,j] <- orth_vec(w[,j], cbind(w[,1:(j-1)]))
    }
    w[,j] <- w[,j]/norm_vec(w[,j])
  }
#############################################
w %<>% abs
w %<>% as.data.frame()
w %<>% mutate(SNP=rownames(w))
w <- merge(w, info, by='SNP')

maf_plot <- function(PC=1){
  comp=paste0('PC', PC)
  arr_t <- desc(w[,comp])
  dat = arrange(w, arr_t)[1:100,]
  score = sum(dat[,comp]*dat$MAF)
  snp_t <- reorder(dat$SNP, desc(dat[,comp]))
  ggplot(data = dat, aes_string(x=snp_t, y=comp)) +
    geom_bar(stat = "identity") + 
    geom_line(aes(x=SNP, y=MAF, group=1))+
    annotate("text", x=40, y=0.15, label= paste('maf score=',round(score, digits = 2)), size=5)
}

p <- lapply(1:15, maf_plot)

do.call("grid.arrange", c(p, ncol=5))

