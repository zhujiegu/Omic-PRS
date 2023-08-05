#!/bin/sh

#SBATCH --time=8:00:00
#SBATCH --mem=10g
#SBATCH --cpus-per-task=1
#SBATCH --output=./rx_mtb/R-%x_screen.out

ulimit -u 20000

rx=$1

# try different nr of rx
#Rscript PO2PLS_gly_fit_rx.R $rx
Rscript PO2PLS_mtb_fit_rx.R $rx

echo 'All complete'

