library(data.table)
library(profvis)
######## SNP and Gene mapping #####
library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb) # was eerst library(GenomeInfo)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(magrittr)
library(org.Hs.eg.db)
library(pryr)
######## Nr of cores ####
library(parallel)
Nr_core=20

for(k in 1:23){
  # gather each chr
  dat<- fread(file = paste('chr',k,'.raw', sep = ''), header = T, verbose = F, nThread = Nr_core)
  
  ## --------------------------------------------------------------------------------
  # remove meta data from matrix
  sample_id <- dat$FID
  dat <- dat[,-(1:6)]
  
  # get position from colname
  snp_names <- colnames(dat)
  snp_pos <- sapply(snp_names, function(e) strsplit(e, '_')[[1]][1])
  
  # check duplicate and remove
  rrm <- which(duplicated(snp_pos))
  snp_pos <- snp_pos[-rrm]
  snp_names <- snp_names[-rrm]
  dat <- dat[ ,-rrm, with = F]
  
  
  ## --------------------------------------------------------------------------------
  
  string2range <- function(pos, delim=' ', region=TRUE) {
    posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
    posp[,1] <- posp[,1]
    posp[,2] <- as.numeric(as.character(posp[,2]))
    if(region) {
      posp[,3] <- as.numeric(as.character(posp[,3]))
    } else {
      posp[,3] <- posp[,2]
    }
    return(posp)
  }
  
  range2GRanges <- function(df) {
    require(GenomicRanges)
    require(IRanges)
    gr <- GenomicRanges::GRanges(
      seqnames = paste('chr',df[,1], sep = ''),
      ranges=IRanges(start = df[,2], end = df[,3])
    )
    return(gr)
  }
  
  snps_ranges <- string2range(snp_pos, delim=":", region=FALSE)
  # head(snps_ranges)
  
  snps_granges <- range2GRanges(snps_ranges)
  names(snps_granges) <- snp_pos
  # head(snps_granges)
  
  
  
  ## ----Match SNPs Genes, warning=FALSE---------------------------------------------
  genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  
  # genes
  
  ##find the "nearest" genes
  genes2 <- genes
  extra_width1 <- 50000L
  extra_width2 <- 50000L
  genes2@ranges@start <- genes@ranges@start - extra_width1
  genes2@ranges@width <- genes@ranges@width + extra_width1 + extra_width2
  
  # # hitted genes / all genes
  # sum(countOverlaps(genes2, snps_granges)>0) / length(genes2)
  
  # matched SNPs / all SNPs
  print(sum(countOverlaps(snps_granges, genes2)>0))
  print(length(snps_granges))
  print(sum(countOverlaps(snps_granges, genes2)>0) / length(snps_granges))
  
  
  ## ----cbind SNPs Genes------------------------------------------------------------
  overlapping_indices <- as.matrix(findOverlaps(snps_granges, genes2))
  geneList <- genes[overlapping_indices[,2]]
  geneList
  
  ##Annotate
  #select(org.Hs.eg.db, mcols(geneList)$gene_id, keytype="ENTREZID", colum=c("SYMBOL"))
  
  ## Cbind with rs values
  rsToGene <- cbind(RefSNP_id=rownames(mcols(snps_granges[overlapping_indices[,1]])), mcols(geneList))
  uniqueSNPs <- !duplicated(rsToGene$RefSNP_id)
  # Now rsToGene contains unique SNPs, ie took the first "nearest" gene
  rsToGene <- rsToGene[uniqueSNPs,]
  rsToGene
  
  
  ## ----list of SNP names per gene--------------------------------------------------
  # Make list of SNP names per gene
  cl <- makePSOCKcluster(detectCores())
  clusterEvalQ(cl, {library(data.table)})
  clusterExport(cl, "rsToGene")
  # This try does not work the first time...
  try(
    system.time(
      gene_list <- parSapply(cl, unique(rsToGene$gene_id), 
                             function(e) {rsToGene[which(rsToGene$gene_id %in% e),][,1]}, USE.NAMES=T)
    ),silent = T)
  
  # so we try again :)
  system.time(
    gene_list <- parSapply(cl, unique(rsToGene$gene_id), 
                           function(e) {rsToGene[which(rsToGene$gene_id %in% e),][,1]}, USE.NAMES=T)
  )
  stopCluster(cl)
  
  
  # ## ----Plot SNPs/Genes per chr, fig.width=14, fig.height=10, warning=FALSE---------
  # ## Plot number of SNPs and genes per chromosome
  # cat("Note:",extra_width1/1000,"Kb added at beginning,",extra_width2/1000, "Kb at end")
  # # take unique SNPs and corresponding genes
  # snpInfo <- cbind(chrom=as.factor(snps_granges[overlapping_indices[,1]][uniqueSNPs]@seqnames), as.data.frame(rsToGene))
  # # make table of number of SNPs on each chromosome, order by chr name
  # tableSNPonChr <- table(snpInfo$chr)[order(as.numeric(stringr::str_sub(names(table(snpInfo$chr)), start=4)))]
  # tableGeneonChr <- table(snpInfo$chr[!duplicated(rsToGene$gene_id)])[order(as.numeric(stringr::str_sub(names(table(snpInfo$chr[!duplicated(rsToGene$gene_id)])), start=4)))]
  # # plot bars with height equal to nr SNPs/genes per chr, first chr is reference height
  # plot(tableSNPonChr, type='h', lwd=1, col='grey',lty=2, ylab='Nr of SNPs/Genes')
  # points(1:23+0.1, tableGeneonChr*tableSNPonChr[1]/tableGeneonChr[1], type='h', col='red', lty=2)
  # text(tableSNPonChr, labels = tableSNPonChr, srt=90, pos = 2)
  # text(1:23+0.1, tableGeneonChr*tableSNPonChr[1]/tableGeneonChr[1], labels = tableGeneonChr, srt=-90, pos = 4, col=2)
  # title(main=paste(sum(tableSNPonChr),"SNPs and",sum(tableGeneonChr), "Genes available and mapped\n",
  #                  "Number of SNPs lost:", sum(!uniqueSNPs)))
  
  
  
  ## ----Plot Genes vs SNPs, fig.width=14, fig.height=10-----------------------------
  ## Plot distribution of number of SNPs per Gene
  cairo_pdf(file = paste("GenevsSNP_chr",k,".pdf",sep = ''), width = 11, height =6, onefile = T, fallback_resolution = 600)
  plot(log10(table(sapply(gene_list, length))), yaxt='n', type='h', lwd=1, ylab='Nr of Genes', xlab='Nr of SNPs on a Gene')
  axis(side=2, at = axTicks(2), labels = signif(10^axTicks(2),4))
  title(main = "Nr of Genes containing certain amount of SNPs")
  dev.off()
  
  
  ## ----Eigenvalus, cache=TRUE------------------------------------------------------
  mean_imput <- function(X){
    for(i in 1:ncol(X)){
      X[which(is.na(X[,i])),i] <- mean(X[,i], na.rm=T)
    }
    return(X)
  }
  
  ## calculate all eigenvalues per gene, takes some time
  system.time(
    eigenvals <- sapply(1:length(gene_list), function(i){
      if(length(which(snp_pos %in% gene_list[[i]]))>0){
        d=svd(scale(mean_imput(as.matrix(dat[,which(snp_pos %in% gene_list[[i]]),with=F]))),0,0)$d
        return(d^2) # / geneList[names(gene_list)[i]]@ranges@width)
      }
      return(numeric(0))
    })
  )
  names(eigenvals) <- names(gene_list)
  SNPsOnGenes <- sapply(eigenvals, length, USE.NAMES = T)
  
  
  
  ## ----Plot Genes/PCs, fig.width=14, fig.height=10---------------------------------
  if(length(which(sapply(eigenvals,length)==0)) !=0 ){
    eigenvals2 <- eigenvals[-(which(sapply(eigenvals,length)==0))]
  }else{
    eigenvals2 <- eigenvals
  }
  
  
  ## plot number of PCs to retain for 80%, 95% variance, then
  #     Nr of SNPs on gene vs Nr of PCs to retain
  PCs_80 <- sapply(eigenvals2, function(e) which(cumsum(e)/sum(e) > 0.8)[1])
  PCs_95 <- sapply(eigenvals2, function(e) which(cumsum(e)/sum(e) > 0.95)[1])
  SNPsPerGene2 <- sapply(eigenvals2, length)
  cairo_pdf(file = paste("NrPC_chr",k,".pdf",sep = ''), width = 11, height =9, onefile = T, fallback_resolution = 600)
  par(mfrow = c(2,2))
  plot(PCs_80 %>% table, ylab=NA, xlab=NA)
  title(main = 'for 80% variance explained', ylab = 'Nr of genes', xlab="Nr of PCs")
  plot(PCs_95 %>% table, ylab=NA, xlab=NA)
  title(main = 'for 95% variance explained', ylab = 'Nr of genes', xlab="Nr of PCs")
  plot(SNPsPerGene2, PCs_80, xlab=NA,ylab=NA)
  title(ylab = 'Nr of PCs', xlab="Nr of SNPs")
  mtext(paste("PCs =", round(coef(lm(PCs_80 ~ SNPsPerGene2))[1],2), "+", round(coef(lm(PCs_80 ~ SNPsPerGene2))[2],3),"* SNPs"), 3, -1.5)
  plot(SNPsPerGene2, PCs_95, xlab=NA,ylab=NA)
  title(ylab = 'Nr of PCs', xlab="Nr of SNPs")
  mtext(paste("PCs =", round(coef(lm(PCs_95 ~ SNPsPerGene2))[1],2), "+", round(coef(lm(PCs_95 ~ SNPsPerGene2))[2],3),"* SNPs"), 3, -1.5)
  par(mfrow = c(1,1))
  dev.off()
  paste("Number of PCs in 80% case:", sum(PCs_80))
  paste("Number of PCs in 95% case:", sum(PCs_95))
  
  
  ## ----SNP-Gene-Score--------------------------------------------------------------
  if(length(which(sapply(eigenvals,length)==0)) !=0 ){
    gene_list2 <- gene_list[-(which(sapply(eigenvals,length)==0))]
  }else{
    gene_list2 <- gene_list
  }
  
  
  print(all.equal(names(PCs_80), names(gene_list2))) # little sanity check
  get_PC1 <- function(i){
    snp_names = gene_list2[[i]]
    snps_on_gene <- which(snp_pos %in% snp_names)
    if(length(snps_on_gene)>0){
      dat_on_gene <- as.matrix(dat[,snps_on_gene,with=F])
      for(j in 1:length(snps_on_gene)){
        dat_on_gene[which(is.na(dat_on_gene[,j])),j] = mean(dat_on_gene[,j], na.rm=T)
      }
      PC1 = svd(scale(dat_on_gene,scale=F), nu=PCs_80[i], nv=0)
      PC1 = PC1$u %*% diag(PC1$d[1:PCs_80[i]], PCs_80[i])
      PC1 = PC1 / PCs_80[i]
      colnames(PC1) <- paste(names(gene_list2)[i],1:PCs_80[i],sep="_")
    } else {
      PC1 = NA*as.matrix(1:nrow(dat))
    }
    return(PC1)
  }
  
  get_PCweights <- function(i){
    snp_names = gene_list2[[i]]
    snps_on_gene <- which(snp_pos %in% snp_names)
    if(length(snps_on_gene)>0){
      dat_on_gene <- as.matrix(dat[,snps_on_gene,with=F])
      for(j in 1:length(snps_on_gene)){
        dat_on_gene[which(is.na(dat_on_gene[,j])),j] = mean(dat_on_gene[,j], na.rm=T)
      }
      PC1 = svd(scale(dat_on_gene,scale=F), nv=PCs_80[i], nu=0)
      PC1 = PC1$v
      colnames(PC1) <- paste(names(gene_list2)[i],1:PCs_80[i],sep="_")
      row.names(PC1) <- as.character(snp_pos[snps_on_gene])
    } else {
      PC1 = NA*as.matrix(1:ncol(dat))
    }
    return(PC1)
  }
  
  system.time(
    SNPweights <- mclapply(1:length(PCs_80), get_PCweights, mc.cores=Nr_core)
  )
  names(SNPweights) <- names(gene_list2)
  
  system.time(
    GeneData <- do.call(cbind, sapply(1:length(PCs_80), get_PC1))
  )
  row.names(GeneData) <- sample_id
  
  saveRDS(gene_list2, file = paste("genelist_chr",k,".rds", sep = ''))
  saveRDS(GeneData, file = paste("GeneData_chr",k,".rds", sep = ''))
}



