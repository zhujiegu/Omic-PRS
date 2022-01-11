# CWP
chronic widespread musculoskeletal pain analysis with genetics and glycomics

plink: QC from the imputed genetic data to dosage data matrix
check_sample_ID.R: check if samples from all chr match after QC.
SNPtoGene.R: summarizes a dosage data matrix to genes (contained in Togene_all.R)
Togene_all.R: read dosage data matrix (per chr) and summarize per gene with PCA, keep PCs that explain 80% of variation. Parallel computing on HPC.
