#!/bin/sh

#SBATCH --job-name="datX"
#SBATCH --time=10:00:00
#SBATCH --mem=100g
#SBATCH --cpus-per-task=5
#SBATCH --array=1-50
#SBATCH --output=./dat_rep/R-%x.%j.out

ulimit -u 20000

# combine chrs

# split samples
#Rscript Training_split.R

# split data to training and test, SNP -> gene
#Rscript SNP_PCA.R
echo 'SNP to gene done'

# keep training, preparing for .ped file
plink2 --bfile chr_all --make-bed --keep ./dat_rep/id_training_${SLURM_ARRAY_TASK_ID}.txt --extract snp_list --out ./dat_rep/chr_all_training_${SLURM_ARRAY_TASK_ID}

# .ped file in training for GWAS
plink --bfile ./dat_rep/chr_all_training_${SLURM_ARRAY_TASK_ID} --recode --tab --out ./dat_rep/chr_all_training_${SLURM_ARRAY_TASK_ID}
echo 'PED file convert complete'

# convert ped file to ABEL raw file
Rscript ABELraw.R
echo 'finish'



