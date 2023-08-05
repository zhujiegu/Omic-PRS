#!/bin/sh

#SBATCH --time=00:30:00
#SBATCH --mem=10g
#SBATCH --cpus-per-task=10
#SBATCH --output=./fit_files_mtb/R-%x.out

ulimit -u 20000

i=$1

# try different nr of rx
Rscript PRSice.R --dir . \
    --prsice PRSice_linux \
    --base /hpc/julius_bs/bs_omics/NMR_GWAS/gwas_${i} --snp ID --chr chromosome --bp position --A1 EA --A2 NEA --stat beta --pvalue p-value \
    --beta \
    --binary-target F \
    --bar-levels 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,0.01,0.1,1 --fastscore \
    --target chr_all_QC_25 \
    --pheno bmi.txt \
    --thread 80 \
    --ignore-fid \
    --all-score \
    --out ./fit_files_mtb/PRS_${i}

echo 'All complete'

