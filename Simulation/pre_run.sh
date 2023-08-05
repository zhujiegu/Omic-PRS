#!/bin/sh

#SBATCH --job-name="pre_run"
#SBATCH --time=10:00:00
#SBATCH --mem=100g
#SBATCH --cpus-per-task=5
#SBATCH --output=R-pre_run.log

ulimit -u 20000

# SNP to genes for all samples
Rscript SNP_PCA_allsample.R

# PCA of SNP and gene data, get loadings as base for true joint loading
Rscript PCA_x.R

# Cut loadings at a percentile to create different causal proportion
Rscript cut_loading.R
