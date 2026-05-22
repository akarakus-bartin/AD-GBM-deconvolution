# =============================================================================
# 01_data_acquisition.R
# Phase 1: Acquisition and preprocessing of bulk transcriptomic datasets
# =============================================================================
# 
# Datasets:
#   - GSE48350 (Berchtold et al.): AD hippocampus + entorhinal cortex,
#                                  Affymetrix HG-U133 Plus 2.0, 75 samples
#   - GSE36980 (Hokama et al.):    AD hippocampus + temporal cortex,
#                                  Affymetrix HuGene 1.0 ST, 47 samples
#   - TCGA-GBM + GTEx (UCSC Xena TOIL): GBM primary tumors + normal brain,
#                                       RNA-seq harmonized, 339 samples
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(GEOquery)
  library(affy)
  library(oligo)
  library(hgu133plus2.db)
  library(hugene10sttranscriptcluster.db)
  library(tidyverse)
})

set.seed(42)

# --- Output directory structure ---
data_raw_dir <- here("data", "raw")
data_proc_dir <- here("data", "processed")
if (!dir.exists(data_raw_dir))  dir.create(data_raw_dir, recursive = TRUE)
if (!dir.exists(data_proc_dir)) dir.create(data_proc_dir, recursive = TRUE)

# =============================================================================
# 1. GSE48350 — Affymetrix HG-U133 Plus 2.0 (Berchtold et al. 2008)
# =============================================================================
cat("=== Downloading GSE48350 (CEL files via GEOquery) ===\n")

gse48350_dir <- file.path(data_raw_dir, "GSE48350")
if (!dir.exists(gse48350_dir)) dir.create(gse48350_dir)

# Download series matrix for metadata
gse48350 <- getGEO("GSE48350", destdir = gse48350_dir,
                   GSEMatrix = TRUE, getGPL = FALSE)
gse48350 <- gse48350[[1]]
meta_48350 <- pData(gse48350)

# Filter to hippocampus and entorhinal cortex only
relevant_regions <- c("hippocampus", "entorhinal cortex")
meta_48350 <- meta_48350[grepl(paste(relevant_regions, collapse = "|"),
                               meta_48350$`brain region:ch1`,
                               ignore.case = TRUE), ]

# Download raw CEL files for selected samples
getGEOSuppFiles("GSE48350", baseDir = data_raw_dir,
                makeDirectory = FALSE, fetch_files = TRUE)
untar(file.path(data_raw_dir, "GSE48350_RAW.tar"), exdir = gse48350_dir)

# RMA normalization
cel_files <- list.files(gse48350_dir, pattern = "\\.CEL\\.gz$",
                        full.names = TRUE)
affy_data <- ReadAffy(filenames = cel_files)
eset_48350 <- affy::rma(affy_data)
expr_48350 <- exprs(eset_48350)

# Probe-to-gene mapping (hgu133plus2.db)
probe_ids <- rownames(expr_48350)
probe_to_symbol <- AnnotationDbi::select(hgu133plus2.db,
                                         keys = probe_ids,
                                         columns = "SYMBOL",
                                         keytype = "PROBEID")
probe_to_symbol <- probe_to_symbol[!is.na(probe_to_symbol$SYMBOL), ]

# Aggregate probes to gene level — keep probe with highest mean expression
probe_means <- rowMeans(expr_48350)
probe_to_symbol$mean_expr <- probe_means[match(probe_to_symbol$PROBEID,
                                               names(probe_means))]
probe_to_symbol <- probe_to_symbol %>%
  group_by(SYMBOL) %>%
  slice_max(mean_expr, n = 1) %>%
  ungroup()

expr_48350_gene <- expr_48350[probe_to_symbol$PROBEID, ]
rownames(expr_48350_gene) <- probe_to_symbol$SYMBOL

# Low-expression gene filter (expressed > 4 in >= 20% of samples)
keep_genes <- rowSums(expr_48350_gene > 4) >= (0.2 * ncol(expr_48350_gene))
expr_48350_gene <- expr_48350_gene[keep_genes, ]

cat("GSE48350 dimensions:", dim(expr_48350_gene), "\n")

# Save
saveRDS(list(expression = expr_48350_gene, metadata = meta_48350),
        file.path(data_proc_dir, "GSE48350_RMA_normalized.rds"))

# =============================================================================
# 2. GSE36980 — Affymetrix HuGene 1.0 ST (Hokama et al. 2014)
# =============================================================================
cat("\n=== Downloading GSE36980 (oligo package for HuGene 1.0 ST) ===\n")

gse36980_dir <- file.path(data_raw_dir, "GSE36980")
if (!dir.exists(gse36980_dir)) dir.create(gse36980_dir)

gse36980 <- getGEO("GSE36980", destdir = gse36980_dir,
                   GSEMatrix = TRUE, getGPL = FALSE)
gse36980 <- gse36980[[1]]
meta_36980 <- pData(gse36980)

getGEOSuppFiles("GSE36980", baseDir = data_raw_dir,
                makeDirectory = FALSE)
untar(file.path(data_raw_dir, "GSE36980_RAW.tar"), exdir = gse36980_dir)

cel_files_36980 <- list.files(gse36980_dir, pattern = "\\.CEL\\.gz$",
                              full.names = TRUE)
oligo_data <- oligo::read.celfiles(cel_files_36980)
eset_36980 <- oligo::rma(oligo_data)
expr_36980 <- exprs(eset_36980)

# Probe-to-gene mapping (hugene10sttranscriptcluster.db)
probe_ids_36980 <- rownames(expr_36980)
probe_map_36980 <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                         keys = probe_ids_36980,
                                         columns = "SYMBOL",
                                         keytype = "PROBEID")
probe_map_36980 <- probe_map_36980[!is.na(probe_map_36980$SYMBOL), ]
probe_means_36980 <- rowMeans(expr_36980)
probe_map_36980$mean_expr <- probe_means_36980[match(probe_map_36980$PROBEID,
                                                     names(probe_means_36980))]
probe_map_36980 <- probe_map_36980 %>%
  group_by(SYMBOL) %>%
  slice_max(mean_expr, n = 1) %>%
  ungroup()

expr_36980_gene <- expr_36980[probe_map_36980$PROBEID, ]
rownames(expr_36980_gene) <- probe_map_36980$SYMBOL

keep_36980 <- rowSums(expr_36980_gene > 4) >= (0.2 * ncol(expr_36980_gene))
expr_36980_gene <- expr_36980_gene[keep_36980, ]

cat("GSE36980 dimensions:", dim(expr_36980_gene), "\n")

saveRDS(list(expression = expr_36980_gene, metadata = meta_36980),
        file.path(data_proc_dir, "GSE36980_RMA_normalized.rds"))

# =============================================================================
# 3. TCGA-GBM + GTEx — via UCSC Xena TOIL pipeline
# =============================================================================
# Note: Pre-processed log2(TPM+0.001) values downloaded from:
#   https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
# This pipeline uniformly aligns and quantifies TCGA + GTEx RNA-seq data,
# removing batch effects related to different processing pipelines.
# =============================================================================
cat("\n=== Loading TCGA-GBM + GTEx (UCSC Xena TOIL) ===\n")

# Expression matrix and phenotype file should be downloaded manually
# from UCSC Xena and placed in: data/raw/xena/
xena_dir <- file.path(data_raw_dir, "xena")
if (!dir.exists(xena_dir)) dir.create(xena_dir)

xena_expr_file <- file.path(xena_dir, "TcgaTargetGtex_RSEM_Hugo_norm_count.gz")
xena_pheno_file <- file.path(xena_dir, "TcgaTargetGTEX_phenotype.txt.gz")

if (file.exists(xena_expr_file) && file.exists(xena_pheno_file)) {
  
  # Load full Xena matrix and filter to GBM + brain GTEx samples
  expr_xena <- data.table::fread(xena_expr_file, data.table = FALSE)
  rownames(expr_xena) <- expr_xena[, 1]
  expr_xena <- expr_xena[, -1]
  
  pheno_xena <- data.table::fread(xena_pheno_file, data.table = FALSE)
  
  # Filter: primary GBM tumors + GTEx hippocampus and frontal cortex (BA9)
  gbm_samples <- pheno_xena %>%
    filter(detailed_category == "Glioblastoma Multiforme",
           `_sample_type` == "Primary Tumor")
  
  normal_samples <- pheno_xena %>%
    filter(`_study` == "GTEX",
           grepl("Hippocampus|Brain - Frontal Cortex \\(BA9\\)",
                 `primary disease or tissue`, ignore.case = TRUE))
  
  selected_samples <- c(gbm_samples$sample, normal_samples$sample)
  selected_samples <- intersect(selected_samples, colnames(expr_xena))
  
  expr_gbm <- expr_xena[, selected_samples]
  pheno_gbm <- pheno_xena[match(selected_samples, pheno_xena$sample), ]
  
  # Add diagnosis label and brain region
  pheno_gbm$diagnosis <- ifelse(pheno_gbm$`_sample_type` == "Primary Tumor",
                                "GBM", "Control")
  pheno_gbm$brain_region <- ifelse(pheno_gbm$diagnosis == "GBM", "tumor",
                                   ifelse(grepl("Hippocampus",
                                                pheno_gbm$`primary disease or tissue`),
                                          "hippocampus", "frontal_cortex"))
  pheno_gbm$sex <- tolower(pheno_gbm$`_gender`)
  
  # Filter low-expression genes
  keep_gbm <- rowSums(expr_gbm > 1) >= (0.2 * ncol(expr_gbm))
  expr_gbm <- expr_gbm[keep_gbm, ]
  
  cat("TCGA-GBM + GTEx dimensions:", dim(expr_gbm), "\n")
  cat("Diagnosis distribution:\n")
  print(table(pheno_gbm$diagnosis))
  
  saveRDS(list(expression = expr_gbm, metadata = pheno_gbm),
          file.path(data_proc_dir, "TCGA_GTEx_GBM_filtered.rds"))
  
} else {
  cat("INFO: UCSC Xena files not found. Download manually from:\n")
  cat("  https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx\n")
  cat("  Files needed:\n")
  cat("    - TcgaTargetGtex_RSEM_Hugo_norm_count.gz\n")
  cat("    - TcgaTargetGTEX_phenotype.txt.gz\n")
  cat("  Place in:", xena_dir, "\n")
}

cat("\n=== Phase 1: Data acquisition complete ===\n")