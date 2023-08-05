# For details of the methods, simulation, and data analysis, please visit Chapter 3 of my PhD thesis: https://dspace.library.uu.nl/handle/1874/428052

## Simulation
### R files:
- `Training_split.R`: split training and test
-	`cut_loading.R`: hard thresholding the eigen vector of genetic data
-	`maf.R`: plot MAF
-	`PCA_x.R`: perform PCA on genetic data and take the first 20 PCs
-	`SNP_PCA.R`: summarise SNP data to GPC data
-	`training_part(_snp).R`: train O2PLS models, perform GWAS for omic-PRS on GPC data (on SNP data)
-	`test_part(_snp).R`: compute metrics to evaluate performance of all methods (on SNP data)
### sh files for executing on HPC (slurm)
-	`pre_run.sh`: run SNP_PCA.R, PCA_x.R, cut_loading.R
-	`simu_datX.sh`: run Training_split.R, SNP_PCA.R, prepare .ped file for GWAS
-	`simu_run(_snp).sh`: Iterate over simulation settings, in each setting, run training_part.R, compute PRSs, and run test_part.R

## Data analysis
### `Analysis_steps`: steps for analysis, containing some plink and PRSice commands
### R files:
- `gly_process.R`, `mtb_process.R`: read raw glycomics and metabolomics data and preprocess
-	`bmi_explore.R`, `db2_explore.R`: preprocessing the outcome BMI, Type II diabetes
- `GWAS_stats_process.R`, `GWAS_stats_process_bash`: formatting downloaded GWAS
-	`overlapping_id.R`: extract sample with omics and outcome, split into training and test
-	`combine_gene_gly`: preprocess genetic and glycan data, ready for O2PLS
-	`SNP_PCA.R`: summarise SNP data to GPC data
-	`gly_bmi.R`: select glycan peaks that are significantly associated with BMI
-	`Nr_comp_bmi.R`: determine and number of components for O2PLS methods
-	`bmi_predict.R`: performance of predicting BMI of all the methods
-	`O2PLS_gly_fit_1000.R`: fit O2PLS methods (without PO2PLS)
-	`PO2PLS_gly_fit_1000.R`: fit PO2PLS
-	`PO2PLS_fit_rx.R`: investigating the impact of the number of x-specific component
### sh files:
-	`run_Nr_comp.sh`: run Nr_comp_bmi.R
-	`run_O2PLS.sh`, `run_PO2PLS.sh`: run O2PLS, PO2PLS.
-	`run_rx.sh`, `rx_sbatch.sh`: investigate the impact of the number of x-specific component
-	`run_PRS.sh`, `sbatch_PRS.sh`: compute PRSs for each glycan
