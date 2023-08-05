library(OmicsPLS)

w_g <- readRDS('PCA_gene.rds')
w_g <- w_g$rotation[,1:20]

w_g_l <- vector(mode = 'list', length = 4)
# w_g_l_tt <- vector(mode = 'list', length = 4)
# 
# for(cut in 1:4){
#   w <- w_g
#   cut_pct <- c(1,0.1,0.01,0.001)[cut]
#   for(j in 1:20){
#     w[abs(w[,j])<quantile(abs(w[,j]),(1-cut_pct)),j] <- 0
#     if(j%%5!=1){
#       w[,j] <- orth_vec(w[,j], cbind(w[,(5*((j-1)%/%5)+1):(j-1)]))
#     }
#     w[,j] <- w[,j]/norm_vec(w[,j])
#   }
#   w_g_l_tt[[cut]] <- w
# }
# 
# crossprod(w_g_l[[4]], w_g_l_tt[[4]])

for(cut in 1:4){
  w <- w_g
  cut_pct <- c(1,0.1,0.01,0.001)[cut]
  for(j in 1:20){
    w[abs(w[,j])<quantile(abs(w[,j]),(1-cut_pct)),j] <- 0
    if(j>1){
      w[,j] <- orth_vec(w[,j], cbind(w[,1:(j-1)]))
    }
    w[,j] <- w[,j]/norm_vec(w[,j])
  }
  w_g_l[[cut]] <- w
}




names(w_g_l) <- c(1,0.1,0.01,0.001)

crossprod(w_g_l[[2]], w_g) %>% diag
# crossprod(w_g_l_tt[[2]], w_g) %>% diag

saveRDS(w_g_l, file = 'w_cut_gene.rds')

non_z <- lapply(1:20, function(e) rownames(w_g)[which(w_g_l[[2]][,e]!=0)])
#############

w_s <- readRDS('PCA_snp.rds')
w_s <- w_s$rotation[,1:20]

w_s_l <- vector(mode = 'list', length = 4)
for(cut in 1:4){
  w <- w_s
  cut_pct <- c(1,0.1,0.01,0.001)[cut]
  if(cut_pct != 0.001){
    for(j in 1:20){
      w[abs(w[,j])<quantile(abs(w[,j]),(1-cut_pct)),j] <- 0
      if(j>1){
        w[,j] <- orth_vec(w[,j], cbind(w[,1:(j-1)]))
      }
      w[,j] <- w[,j]/norm_vec(w[,j])
    }
  }else{  # a bug in OmicsPLS, fix here temporarily
    for(j in 1:20){
      w[abs(w[,j])<quantile(abs(w[,j]),(1-cut_pct)),j] <- 0
      if(j>1&j!=15){
        w[,j] <- orth_vec(w[,j], cbind(w[,1:(j-1)]))
      }
      w[,j] <- w[,j]/norm_vec(w[,j])
    }
  }

  w_s_l[[cut]] <- w
}
names(w_s_l) <- c(1,0.1,0.01,0.001)
saveRDS(w_s_l, file = 'w_cut_snp.rds')

# 
# crossprod(w_s_l[[4]], w_s) %>% diag
# crossprod(w_s_l[[4]])

non_z <- lapply(1:20, function(e) rownames(w_s)[which(w_s_l[[3]][,e]!=0)])

print('loading cutting done')