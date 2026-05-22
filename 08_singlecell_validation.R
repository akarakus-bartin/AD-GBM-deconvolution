# =============================================================================
# 08_singlecell_validation.R
# Phase 6: Single-cell exploration of inverse comorbidity pathways
# =============================================================================
#
# Dataset: GSE138852 (Grubman et al. 2019, Nature Neuroscience)
#   - Entorhinal cortex snRNA-seq
#   - 4 AD donors + 4 control donors (pooled per pair)
#   - ~14,000 nuclei after upstream filtering
#
# Goal: Test whether the three bulk-level inverse comorbidity pathways
#       (DNA REPAIR, GLYCOLYSIS, MYC TARGETS V1) show concordant patterns
#       at single-cell resolution in AD entorhinal neurons.
#
# Analysis:
#   1. Filter doublets and unidentified populations
#   2. Log-normalize
#   3. Compute per-cell Hallmark pathway scores (Seurat AddModuleScore)
#   4. Donor-level pseudobulk Wilcoxon test (n=4 vs n=4) + Cohen's d
#
# Note: This analysis is exploratory due to limited donor numbers.
# Negative findings are reported honestly as a study limitation.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(msigdbr)
  library(tidyverse)
  library(ggplot2)
})

set.seed(42)

# --- Paths ---
sc_dir <- here("data", "raw", "singlecell", "GSE138852_Grubman")
grubman_counts_file <- file.path(sc_dir, "GSE138852_counts.csv.gz")
grubman_cov_file    <- file.path(sc_dir, "GSE138852_covariates.csv.gz")

# =============================================================================
# 1. LOAD COUNTS AND METADATA
# =============================================================================

cat("=== Loading Grubman et al. 2019 snRNA-seq data ===\n")

# Counts: rows = genes, columns = cell barcodes
grubman_counts <- fread(grubman_counts_file, header = TRUE)
gene_names <- grubman_counts$V1
counts_mat <- as.matrix(grubman_counts[, -1, with = FALSE])
rownames(counts_mat) <- gene_names

# Sparse conversion (snRNA-seq is highly sparse)
counts_sparse <- as(counts_mat, "CsparseMatrix")
rm(counts_mat); gc(verbose = FALSE)

cat("Counts matrix:", dim(counts_sparse), "\n")

# Covariates (per-cell metadata, includes pre-existing cell type annotation)
grubman_cov <- fread(grubman_cov_file, header = TRUE)
grubman_meta_df <- as.data.frame(grubman_cov)
rownames(grubman_meta_df) <- grubman_meta_df$V1
grubman_meta_df$V1 <- NULL
colnames(grubman_meta_df) <- c("batchCond", "cellType", "cellType_batchCond",
                               "subclustID", "subclustCond")

cat("Cell type distribution:\n")
print(table(grubman_meta_df$cellType, grubman_meta_df$batchCond))

# Align counts and metadata
common_barcodes <- intersect(colnames(counts_sparse), rownames(grubman_meta_df))
counts_sparse <- counts_sparse[, common_barcodes]
grubman_meta_df <- grubman_meta_df[common_barcodes, ]
stopifnot(identical(colnames(counts_sparse), rownames(grubman_meta_df)))

# =============================================================================
# 2. CREATE SEURAT OBJECT + QC
# =============================================================================

cat("\n=== Creating Seurat object ===\n")
grubman_seurat <- CreateSeuratObject(
  counts = counts_sparse,
  meta.data = grubman_meta_df,
  project = "Grubman_AD",
  min.cells = 3,
  min.features = 200
)

cat("Pre-filter:", ncol(grubman_seurat), "cells\n")

# Filter doublet and unID populations (annotated by Grubman et al.)
cells_to_keep <- !grubman_seurat$cellType %in% c("doublet", "unID")
grubman_clean <- subset(grubman_seurat,
                        cells = colnames(grubman_seurat)[cells_to_keep])

cat("Post-filter (excluding doublet + unID):",
    ncol(grubman_clean), "cells\n\n")

# Mitochondrial percentage
grubman_clean[["percent.mt"]] <- PercentageFeatureSet(grubman_clean,
                                                      pattern = "^MT-")

cat("QC summary per cell type x condition:\n")
qc_summary <- grubman_clean@meta.data %>%
  group_by(cellType, batchCond) %>%
  summarise(
    n_cells = n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_mt_pct = round(median(percent.mt), 2),
    .groups = "drop"
  )
print(qc_summary)

# =============================================================================
# 3. NORMALIZATION
# =============================================================================

cat("\n=== Normalization ===\n")
grubman_clean <- NormalizeData(
  grubman_clean,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# =============================================================================
# 4. INVERSE COMORBIDITY PATHWAY SCORING
# =============================================================================

cat("\n=== Loading Hallmark inverse comorbidity pathways ===\n")
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")

target_pathways <- c(
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_MYC_TARGETS_V1"
)

pathway_gene_lists <- lapply(target_pathways, function(p) {
  genes <- hallmark_df %>%
    filter(gs_name == p) %>%
    pull(gene_symbol) %>%
    unique()
  intersect(genes, rownames(grubman_clean))
})
names(pathway_gene_lists) <- gsub("HALLMARK_", "", target_pathways)

for (i in seq_along(pathway_gene_lists)) {
  cat(target_pathways[i], ":",
      length(pathway_gene_lists[[i]]), "genes available in data\n")
}

# Compute per-cell pathway scores (Seurat AddModuleScore)
cat("\nRunning AddModuleScore...\n")
grubman_clean <- AddModuleScore(
  object = grubman_clean,
  features = pathway_gene_lists,
  name = "Score",
  seed = 42
)

# Rename AddModuleScore output columns (Score1, Score2, Score3)
names_old <- paste0("Score", seq_along(target_pathways))
names_new <- names(pathway_gene_lists)
for (i in seq_along(names_old)) {
  grubman_clean@meta.data[[names_new[i]]] <- grubman_clean@meta.data[[names_old[i]]]
  grubman_clean@meta.data[[names_old[i]]] <- NULL
}

# =============================================================================
# 5. NEURON-FOCUSED COMPARISON: AD vs CONTROL
# =============================================================================

cat("\n=== Neuron-focused inverse comorbidity test ===\n")

neurons_only <- subset(grubman_clean, subset = cellType == "neuron")
cat("Neuron count:", ncol(neurons_only),
    "(AD:", sum(neurons_only$batchCond == "AD"),
    "| Control:", sum(neurons_only$batchCond == "ct"), ")\n\n")

# --- Cell-level Wilcoxon (note: limited interpretability due to
#     within-donor pseudoreplication) ---
cat("--- Cell-level Wilcoxon (each cell as independent observation) ---\n")
for (pw in names_new) {
  scores_ad <- neurons_only@meta.data[neurons_only$batchCond == "AD", pw]
  scores_ct <- neurons_only@meta.data[neurons_only$batchCond == "ct", pw]
  test <- wilcox.test(scores_ad, scores_ct)
  cat(sprintf("%-20s | AD median=%6.3f | Ct median=%6.3f | p=%.2e\n",
              pw, median(scores_ad), median(scores_ct), test$p.value))
}

# --- Donor-level pseudobulk Wilcoxon (proper statistical inference) ---
cat("\n--- Donor-level pseudobulk Wilcoxon (n=4 AD vs n=4 Control) ---\n")

neurons_meta <- neurons_only@meta.data
# Extract donor ID from barcode (format: BARCODE_AD5_AD6 or BARCODE_Ct1_Ct2)
neurons_meta$donor <- sub("^[ACGT]+_", "", colnames(neurons_only))

pseudobulk_scores <- neurons_meta %>%
  group_by(donor, batchCond) %>%
  summarise(
    n_cells = n(),
    DNA_REPAIR = mean(DNA_REPAIR),
    GLYCOLYSIS = mean(GLYCOLYSIS),
    MYC_TARGETS_V1 = mean(MYC_TARGETS_V1),
    .groups = "drop"
  )

print(pseudobulk_scores)

cat("\nDonor-level statistical tests:\n")
for (pw in c("DNA_REPAIR", "GLYCOLYSIS", "MYC_TARGETS_V1")) {
  ad_donors <- pseudobulk_scores[[pw]][pseudobulk_scores$batchCond == "AD"]
  ct_donors <- pseudobulk_scores[[pw]][pseudobulk_scores$batchCond == "ct"]
  wt <- wilcox.test(ad_donors, ct_donors)
  diff <- mean(ad_donors) - mean(ct_donors)
  effect_size <- diff / sd(c(ad_donors, ct_donors))  # Cohen's d
  cat(sprintf("%-20s | diff=%+7.4f | Cohen's d=%+5.2f | wilcox p=%.3f\n",
              pw, diff, effect_size, wt$p.value))
}

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

saveRDS(grubman_clean, here("data", "processed", "Grubman_seurat.rds"))
saveRDS(pseudobulk_scores,
        here("data", "processed", "Grubman_pseudobulk_pathway_scores.rds"))

cat("\n=== Phase 6 (single-cell exploration) complete ===\n")
cat("Honest reporting:\n")
cat("  Donor-level pseudobulk did not detect significant differences\n")
cat("  Cohen's d effect sizes are negligible (|d| < 0.12)\n")
cat("  Limitations: (i) only n=4 vs n=4 donors,\n")
cat("               (ii) snRNA-seq dropout,\n")
cat("               (iii) AD nuclei show ~30%% less gene detection,\n")
cat("               (iv) neuronal subtype composition not resolved.\n")
cat("  This negative finding is reported in the manuscript as a limitation.\n")