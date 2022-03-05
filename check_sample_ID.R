# HPC Plink2 files
sam_id <- list(length=23)
for (i in 1:23) {
  sam_id[i] <- read.table(paste('chr',i,'.psam',sep = ''))
}

for(i in 1:22){
  print(i)
  print(all.equal(sam_id[[i]], sam_id[[i+1]]))
}


# HPC Plink1 files
sam_id <- list(length=22)
for (i in 1:22) {
  sam_id[i] <- read.table(paste('chr',i,'.fam',sep = ''))
}

for(i in 1:21){
  print(i)
  print(all.equal(sam_id[[i]], sam_id[[i+1]]))
}