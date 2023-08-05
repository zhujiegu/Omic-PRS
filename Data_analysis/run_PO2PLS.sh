#!/bin/sh

#SBATCH --job='PO2_screen_scaled_1000'
#SBATCH --time=30:00:00
#SBATCH --mem=20g
#SBATCH --cpus-per-task=1
#SBATCH --output=R-%x.out

ulimit -u 20000

Rscript PO2PLS_gly_fit_1000.R
Rscript PO2PLS_mtb_fit_1000.R

echo 'All complete'

