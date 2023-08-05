#!/bin/sh

for i in `cat /hpc/julius_bs/bs_omics/NMR_GWAS/mtb_sel_names`
do 
	sbatch -J 'PRS_'$i sbatch_PRS.sh $i
done

