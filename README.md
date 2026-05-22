# AD-GBM Comparative Transcriptomics with Cell-Type Deconvolution

Reproducible analysis code for:

> Karakuş A. **Apparent shared signature between Alzheimer's disease and glioblastoma reflects converging neuronal loss, not molecular convergence.** *Submitted to Briefings in Bioinformatics*, 2026.

## Overview

This repository contains R scripts implementing a cell-type-deconvolution-controlled pipeline that systematically re-evaluates apparent transcriptomic convergence between Alzheimer's disease (AD) and glioblastoma multiforme (GBM).

**Central finding:** The 198 shared differentially expressed genes between AD and GBM detected by conventional bulk analysis (96.5% directionally concordant) collapse to **zero** when cell-type composition is controlled. The entirety of apparent AD-GBM transcriptomic convergence reflects shared neuronal loss and glial/immune infiltration rather than convergent cell-autonomous molecular programs.

The pipeline is applied across five complementary analytical layers — gene-level differential expression, co-expression network analysis (WGCNA), module hub gene identification, ranked pathway enrichment (fgsea + Hallmark), and transcription factor activity inference (decoupleR + CollecTRI) — to test whether any cell-composition-resistant signal exists between the two disorders. Three Hallmark pathways (DNA repair, glycolysis, MYC targets) and 14 transcription factors (including OLIG2, MAX, TCF7, HIPK2) are identified that diverge in opposite directions between cell-composition-controlled AD and GBM signatures, providing the first transcriptomic-level evidence consistent with the AD-cancer inverse epidemiological comorbidity.

## Datasets

| Source | Cohort | Platform | n | Phase |
|---|---|---|---|---|
| GSE48350 (GEO) | AD hippocampus + entorhinal cortex | Affymetrix HG-U133 Plus 2.0 | 75 | Bulk AD |
| GSE36980 (GEO) | AD hippocampus + temporal cortex | Affymetrix HuGene 1.0 ST | 47 | Bulk AD |
| TCGA-GBM + GTEx (UCSC Xena TOIL) | GBM primary tumors + normal brain | RNA-seq harmonized | 339 | Bulk GBM |
| GSE138852 (Grubman et al. 2019) | AD entorhinal cortex snRNA-seq | 10x Genomics | 4 AD + 4 Ct | Single-cell |

Meta-cohort assembly combines GSE48350 and GSE36980 (n = 122 total) with ComBat batch correction.

## Pipeline structure

| Script | Phase | Purpose |
|---|---|---|
| `01_data_acquisition.R` | 1 | Download and preprocess raw datasets |
| `02_qc_visualization.R` | 1 QC | Quality control plots (PCA, density, distance) |
| `03_cell_deconvolution.R` | 2 | BRETIGEA cell-type scoring (6 brain cell types) |
| `04_deg_analysis.R` | 3 | limma DEG (uncorrected + deconvolution-controlled) |
| `05_meta_analysis.R` | 3.5 | ComBat batch correction + meta-cohort DEG |
| `06_wgcna_analysis.R` | 4 | Weighted co-expression network analysis |
| `07_gsea_analysis.R` | 5 | Hallmark pathway enrichment (fgsea) |
| `08_singlecell_validation.R` | 6 | Grubman snRNA-seq exploration |
| `09_TF_analysis.R` | 7 | decoupleR + CollecTRI TF activity |
| `10_hub_gene_analysis.R` | 6.5 | WGCNA module hub gene identification |
| `99_publication_figures.R` | Output | Main figures 1–10 |
| `100_supplementary_figures.R` | Output | Supplementary figures S1–S4 |

## Reproducing the analysis

### Software environment

- R 4.6.0 on macOS (Apple Silicon tested)
- Bioconductor 3.23

### Key packages

```r
# CRAN
install.packages(c("here", "tidyverse", "data.table", "ggplot2",
                   "ggrepel", "patchwork", "pheatmap", "RColorBrewer",
                   "car", "Matrix"))

# Bioconductor
BiocManager::install(c("limma", "sva", "fgsea", "msigdbr",
                       "GEOquery", "affy", "oligo",
                       "hgu133plus2.db", "hugene10sttranscriptcluster.db",
                       "WGCNA", "BRETIGEA", "decoupleR", "OmnipathR"))

# Single-cell
install.packages("Seurat")
```

### Workflow

1. Clone this repository
2. Create a project directory structure mirroring `~/AD_GBM_v2/`:
3. Run scripts in numerical order (01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 99 → 100)
4. Each script generates RDS objects that are read by subsequent scripts

## Data availability

Processed RDS objects (normalized expression matrices, sample metadata with BRETIGEA cell-type scores, DEG tables, WGCNA networks with module assignments, hub gene rankings, GSEA results, TF activity matrices) are deposited at Zenodo:

**DOI:** *[To be added after Zenodo upload]*

Raw data sources are publicly available:
- GSE48350, GSE36980, GSE138852: [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/)
- TCGA-GBM + GTEx: [UCSC Xena TOIL](https://toil.xenahubs.net/)

## Key methodological points

- **Cell-type deconvolution** via BRETIGEA (6 cell types: neuron, astrocyte, microglia, oligodendrocyte, OPC, endothelial)
- **Meta-cohort batch correction** via ComBat (PC1×dataset r reduced from −0.998 to −0.06)
- **Multicollinearity assessment** via VIF (diagnosis VIF = 1.46, well below the 5 threshold)
- **Statistical asymmetry analysis**: AD null TF result (0/758 at FDR<0.05) contrasted with GBM strong signal (258/758 at FDR<0.05)
- **Soft inverse threshold** for hypothesis-generating TF identification (relaxed AD p-value criterion clearly distinguished from confirmatory findings)

## Citation

If you use this code or its findings, please cite:

> Karakuş A. (2026). Apparent shared signature between Alzheimer's disease and glioblastoma reflects converging neuronal loss, not molecular convergence. *Briefings in Bioinformatics*, [DOI to be assigned].

## License

This code is released under the [MIT License](LICENSE). You are free to use, modify, and redistribute the code with attribution.

## Contact

Ahmet Karakuş
Bartın University, Bartın, Turkey
ORCID: *[0000-0003-1458-808X]*
Email: *[akarakus@bartin.edu.tr]*

## Acknowledgments

We thank the data depositors of the publicly available datasets used in this study (GSE48350, GSE36980, GSE138852, TCGA-GBM, GTEx).
