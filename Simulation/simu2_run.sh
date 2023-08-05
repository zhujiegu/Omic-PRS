#!/bin/sh

#SBATCH --time=15:00:00
#SBATCH --mem=25g
#SBATCH --cpus-per-task=1
#SBATCH --array=1-20
#SBATCH --output=./results_simu2/R-%x.%j.out

ulimit -u 20000

suffix='simu2_'${SLURM_ARRAY_TASK_ID}
echo $suffix

# generate training and test data, fit O2PLS, perform GWAS
Rscript training_part_simu2.R

# generate PRS for yy based on summary stats in training
{
for i in {1..10}
do
Rscript PRSice.R --dir . \
    --prsice PRSice_linux \
    --base ./fit_files_simu2/GWAS_yy${i}_${suffix}.txt --snp SNP --chr Chromosome --bp Position --A1 A2 --A2 A1 --stat effB --pvalue P2df \
    --beta \
    --lower 1e-10 \
    --bar-levels 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,0.01,0.1,1 --fastscore \
    --target ./dat_rep/chr_all_training_${SLURM_ARRAY_TASK_ID} \
    --thread 1 \
    --binary-target F \
    --pheno ./fit_dat_simu2/yy_prsice_${suffix} --pheno-col yy${i} \
    --ignore-fid \
    --perm 99999 \
    --out ./fit_files_simu2/PRS_yy${i}_${suffix}
done
} &> /dev/null
echo 'PRS YY in training done'

{
Rscript PRSice.R --dir . \
    --prsice PRSice_linux \
    --base ./fit_files_simu2/GWAS_z_${suffix}.txt --snp SNP --chr Chromosome --bp Position --A1 A2 --A2 A1 --stat effB --pvalue P2df \
    --beta \
    --lower 1e-10 \
    --bar-levels 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,0.01,0.1,1 --fastscore \
    --target ./dat_rep/chr_all_training_${SLURM_ARRAY_TASK_ID} \
    --thread 1 \
    --binary-target F \
    --pheno ./fit_dat_simu2/z_prsice_${suffix} \
    --ignore-fid \
    --perm 99999 \
    --out ./fit_files_simu2/PRS_z_${suffix}
} &> /dev/null
echo 'PRS Z in training done'

# generate PRS based on summary stats and best p-value from training in test
{
for i in {1..10}
do
Rscript PRSice.R --dir . \
    --prsice PRSice_linux \
    --base ./fit_files_simu2/GWAS_yy${i}_${suffix}.txt --snp SNP --chr Chromosome --bp Position --A1 A2 --A2 A1 --stat effB --pvalue P2df \
    --beta \
    --no-full --bar-levels $(tail -n 1 ./fit_files_simu2/PRS_yy${i}_${suffix}.summary | awk '{print $12;}') --fastscore \
    --target chr_all \
    --remove ./dat_rep/id_training_${SLURM_ARRAY_TASK_ID}.txt \
    --thread 1 \
    --binary-target F \
    --no-regress \
    --ignore-fid \
    --print-snp \
    --out ./fit_files_simu2/PRS_yy${i}_${suffix}_test
done
} &> /dev/null
echo 'PRS Y in test done'

{
Rscript PRSice.R --dir . \
    --prsice PRSice_linux \
    --base ./fit_files_simu2/GWAS_z_${suffix}.txt --snp SNP --chr Chromosome --bp Position --A1 A2 --A2 A1 --stat effB --pvalue P2df \
    --beta \
    --no-full --bar-levels $(tail -n 1 ./fit_files_simu2/PRS_z_${suffix}.summary | awk '{print $12;}') --fastscore \
    --target chr_all \
    --remove ./dat_rep/id_training_${SLURM_ARRAY_TASK_ID}.txt \
    --thread 1 \
    --binary-target F \
    --no-regress \
    --ignore-fid \
    --print-snp \
    --out ./fit_files_simu2/PRS_z_${suffix}_test
} &> /dev/null
echo 'PRS Z in test done'

# evaluate performance in training and testing
Rscript test_part_simu2.R
echo 'All complete'





