#!/bin/sh

#SBATCH --job='Nr_comp_bmi'
#SBATCH --time=6:00:00
#SBATCH --mem=20g
#SBATCH --cpus-per-task=1
#SBATCH --output=R-%x.out

ulimit -u 20000

#Rscript Nr_comp.R
#Rscript Nr_comp_gly.R
#Rscript Nr_comp_mtb.R
Rscript Nr_comp_bmi.R
