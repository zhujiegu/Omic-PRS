#!/bin/sh

#SBATCH --job='O2_screen5e8_scale'
#SBATCH --time=8:00:00
#SBATCH --mem=100g
#SBATCH --cpus-per-task=1
#SBATCH --output=R-%x.out

ulimit -u 20000

Rscript O2PLS_gly_fit_1000.R
Rscript O2PLS_mtb_fit_1000.R

echo 'All complete'

