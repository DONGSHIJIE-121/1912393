# PANoptosis-associated transcriptomic analysis in sarcopenia

This repository contains the streamlined R workflows used for the central computational analyses in the manuscript. The scripts retain data import, preprocessing, core statistical analysis, feature prioritization, and result export. Trial code, local absolute paths, and figure-color tuning have been removed.

## Repository structure

```text
code/
  01_bulk_preprocessing_DESeq2.R
  02_functional_enrichment.R
  03_machine_learning_and_validation.R
  04_CIBERSORT_immune_infiltration.R
  05_snRNAseq_GSE167186.R
  README_01_bulk_preprocessing_DESeq2.md
  README_02_functional_enrichment.md
  README_03_machine_learning_and_validation.md
  README_04_CIBERSORT_immune_infiltration.md
  README_05_snRNAseq_GSE167186.md
data/input/
  bulk/
  gene_sets/
  immune/
  snrna/
results/
```

## Analysis order

Run all commands from the repository root.

```r
source("code/01_bulk_preprocessing_DESeq2.R")
source("code/02_functional_enrichment.R")
source("code/03_machine_learning_and_validation.R")
source("code/04_CIBERSORT_immune_infiltration.R")
source("code/05_snRNAseq_GSE167186.R")
```

The bulk and single-nucleus workflows are independent. Script 02 and Script 03 use outputs from Script 01. Script 04 uses the training-cohort TPM matrix. Script 05 uses the GSE167186 Seurat object and Scrublet doublet calls.

## Main analytical parameters

- Training datasets: GSE111006, GSE111010, GSE226151, and the human muscle-tissue subset of GSE238215.
- External validation dataset: GSE111016.
- Low-count filter: count >=10 in at least 50% of training samples.
- Batch adjustment: ComBat-seq with dataset as batch and group retained as the biological variable.
- Differential expression: DESeq2 design `~ Dataset + Group`; adjusted P value <0.05 and absolute log2 fold change >0.25.
- LASSO: binomial model, alpha=1, stratified 10-fold cross-validation, lambda.1se.
- Random forest: 300 trees; top 10 variables ranked by MeanDecreaseGini.
- XGBoost: stratified 5-fold cross-validation; early stopping; top 10 variables ranked by gain.
- CIBERSORT: LM22, 1,000 permutations, quantile normalization disabled.
- snRNA-seq: Seurat 5 workflow, patient-specific 3-MAD quality control, Scrublet doublet removal, 2,000 variable genes, first 20 PCs, clustering resolution 0.5, and AddModuleScore using the prespecified top-50 PANoptosis-associated genes.

## Reproducibility notes

All scripts use relative paths and fixed random seeds. Each script writes `sessionInfo.txt` to its result directory. Raw data and licensed CIBERSORT resources are not redistributed. Users must obtain the source data from GEO and the CIBERSORT R script and LM22 matrix from the official CIBERSORT distribution.

The single-nucleus cluster annotation is fixed to the 12-cluster solution reported in the manuscript. If rerunning with a different Seurat version or altered filtering changes cluster numbering, canonical markers must be checked before applying the annotation map.
