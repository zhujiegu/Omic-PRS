# CWP
chronic widespread musculoskeletal pain analysis with genetics and glycomics

R files:
plink: QC from the imputed genetic data to dosage data matrix, or to plink1 format for lassosum
check_sample_ID.R: check if samples from all chr match after QC.
SNPtoGene.R: summarizes a dosage data matrix to genes (contained in Togene_all.R)
Togene_all.R: read dosage data matrix (per chr) and summarize per gene with PCA, keep PCs that explain 80% of variation. Parallel computing on HPC.
combine_chr.R: combine summarized gene data of each chr to one dataset
combine_gene_gly.R: combine summarized gene data with glycan data, save data ready as input for O2PLS
samples_overlap.R: devide samples to 1.samples with only genetics, 2.overlapping genetics and glycomics. Then de-couple twins. output sample ID list.
gly_sel_cwp.R: test which glycans are associated with cwp, and 1. save glycan data containing these glycans. 2. extract summary stats for these glycans for PRSice
covar.R: extract covariates information from quetionaires, inputation, output a dataframe. Age, Sex, Smoking, Alchohol, Bmi included
Nr_comp.R: determine number of components for O2PLS
PRSice2.txt: command to calculate PRS in PRSice2
sample_split_order.R: sample orders of the training and test genetic data in plink2 format (pgen)
lassosum.R: lassosum for PRS for CWP

Data files:
id_decpl: independent sample IDS, one in twin pair removed
id_A_decpl: SNP + glycomics + CWP
id_B_decpl: SNP + CWP
id_C_decpl: SNP
id_AB_cwp: id of id_A_decpl+id_B_decpl and corresponding cwp outcome
