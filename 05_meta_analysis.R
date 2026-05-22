# =============================================================================
# 05_meta_analysis.R
# Phase 3.5: Meta-cohort assembly (GSE48350 + GSE36980) and meta-DEG analysis
# =============================================================================
#
# Steps:
#   1. Merge two AD datasets at gene-symbol level (17,315 common genes)
#   2. ComBat batch correction (sva package) with dataset as batch variable,
#      diagnosis and brain region preserved as biological covariates
#   3. BRETIGEA cell-type deconvolution on the merged matrix
#   4. limma DEG analysis (uncorrected and deconvolution-controlled)
#   5. Cross-cohort intersection with GBM
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(sva)         # ComBat
  library(limma)
  library(BRETIGEA)
  library(tidyverse)
})

set.seed(42)

# --- Load datasets ---
gse48350 <- readRDS(here("data", "processed", "GSE48350_RMA_normalized.rds"))
gse36980 <- readRDS(here("data", "processed", "GSE36980_RMA_normalized.rds"))

expr_48350 <- gse48350$expression
meta_48350 <- gse48350$metadata

expr_36980 <- gse36980$expression
meta_36980 <- gse36980$metadata

# =============================================================================
# 1. METADATA HARMONIZATION
# =============================================================================

# GSE48350 metadata harmonization
meta_48350_clean <- data.frame(
  sample_id    = rownames(meta_48350),
  diagnosis    = ifelse(grepl("Alzheimer", meta_48350$`disease state:ch1`,
                              ignore.case = TRUE), "AD", "Control"),
  brain_region = tolower(gsub("\\s+", "_", meta_48350$`brain region:ch1`)),
  sex          = tolower(meta_48350$`Sex:ch1`),
  age          = as.numeric(meta_48350$`age (yrs):ch1`),
  dataset      = "GSE48350",
  stringsAsFactors = FALSE
)

# GSE36980 metadata harmonization
meta_36980_clean <- data.frame(
  sample_id    = rownames(meta_36980),
  diagnosis    = ifelse(grepl("AD", meta_36980$`source_name_ch1`,
                              ignore.case = TRUE), "AD", "Control"),
  brain_region = tolower(gsub("\\s+", "_",
                              gsub(".*from\\s+", "",
                                   meta_36980$source_name_ch1))),
  sex          = tolower(meta_36980$`Sex:ch1`),
  age          = as.numeric(meta_36980$`age:ch1`),
  dataset      = "GSE36980",
  stringsAsFactors = FALSE
)

# Combine metadata
meta_combined <- rbind(meta_48350_clean, meta_36980_clean)
rownames(meta_combined) <- meta_combined$sample_id

cat("Combined metadata:\n")
print(table(meta_combined$diagnosis, meta_combined$dataset))

# =============================================================================
# 2. EXPRESSION MATRIX MERGE (common genes)
# =============================================================================

common_genes <- intersect(rownames(expr_48350), rownames(expr_36980))
cat("\nCommon genes between datasets:", length(common_genes), "\n")

expr_48350_sub <- expr_48350[common_genes, ]
expr_36980_sub <- expr_36980[common_genes, ]
expr_merged <- cbind(expr_48350_sub, expr_36980_sub)

# Align metadata with merged expression
meta_combined <- meta_combined[colnames(expr_merged), ]
stopifnot(identical(colnames(expr_merged), rownames(meta_combined)))

# =============================================================================
# 3. PRE-COMBAT PCA DIAGNOSTIC
# =============================================================================

cat("\n--- Pre-ComBat PCA diagnostic ---\n")
top_var <- order(apply(expr_merged, 1, var), decreasing = TRUE)[1:2000]
pca_pre <- prcomp(t(expr_merged[top_var, ]), scale. = TRUE)
ve_pre <- pca_pre$sdev^2 / sum(pca_pre$sdev^2) * 100

pc1_dataset_pre <- cor(pca_pre$x[, 1],
                       as.numeric(factor(meta_combined$dataset)))
pc1_diagnosis_pre <- cor(pca_pre$x[, 1],
                         as.numeric(factor(meta_combined$diagnosis)))

cat("PC1 variance:", round(ve_pre[1], 1), "%\n")
cat("PC1 x dataset r:",   round(pc1_dataset_pre, 3),
    "(strong dataset effect expected)\n")
cat("PC1 x diagnosis r:", round(pc1_diagnosis_pre, 3), "\n")

# =============================================================================
# 4. ComBat BATCH CORRECTION
# =============================================================================

cat("\n--- Running ComBat batch correction ---\n")
mod_combat <- model.matrix(~ diagnosis + brain_region, data = meta_combined)

expr_combat <- ComBat(
  dat = as.matrix(expr_merged),
  batch = meta_combined$dataset,
  mod = mod_combat,
  par.prior = TRUE,
  prior.plots = FALSE
)

# Post-ComBat PCA
pca_post <- prcomp(t(expr_combat[top_var, ]), scale. = TRUE)
ve_post <- pca_post$sdev^2 / sum(pca_post$sdev^2) * 100

pc1_dataset_post <- cor(pca_post$x[, 1],
                        as.numeric(factor(meta_combined$dataset)))
pc1_diagnosis_post <- cor(pca_post$x[, 1],
                          as.numeric(factor(meta_combined$diagnosis)))

cat("Post-ComBat PC1 x dataset r:",   round(pc1_dataset_post, 3),
    "(should be near zero)\n")
cat("Post-ComBat PC1 x diagnosis r:", round(pc1_diagnosis_post, 3),
    "(should be stronger)\n")

# =============================================================================
# 5. BRETIGEA CELL-TYPE DECONVOLUTION (on ComBat-corrected matrix)
# =============================================================================

cat("\n--- BRETIGEA deconvolution on meta-cohort ---\n")
cell_scores_meta <- brainCells(
  inputMat = expr_combat,
  nMarker = 50,
  species = "human",
  celltypes = c("ast", "end", "mic", "neu", "oli", "opc"),
  method = "PCA",
  scale = TRUE
)

# Append cell scores to metadata
meta_combined <- cbind(meta_combined, t(cell_scores_meta))

# =============================================================================
# 6. META-COHORT DEG ANALYSIS
# =============================================================================

meta_combined$diagnosis <- factor(meta_combined$diagnosis,
                                  levels = c("Control", "AD"))
meta_combined$brain_region <- factor(meta_combined$brain_region)
meta_combined$sex <- factor(meta_combined$sex)
meta_combined$dataset <- factor(meta_combined$dataset)

# --- Uncorrected model ---
cat("\n--- Meta: Uncorrected limma model ---\n")
design_meta_uncorr <- model.matrix(
  ~ diagnosis + brain_region + sex + age + dataset,
  data = meta_combined
)
fit_uncorr <- lmFit(expr_combat, design_meta_uncorr)
fit_uncorr <- eBayes(fit_uncorr)
meta_uncorr <- topTable(fit_uncorr, coef = "diagnosisAD",
                        number = Inf, adjust.method = "BH", sort.by = "P")
meta_sig_uncorr <- meta_uncorr[meta_uncorr$adj.P.Val < 0.05 &
                                 abs(meta_uncorr$logFC) > 0.5, ]
cat("Meta uncorrected DEGs:", nrow(meta_sig_uncorr), "\n")

# --- Deconvolution-controlled model ---
cat("--- Meta: Deconvolution-controlled limma model ---\n")
design_meta_ctrl <- model.matrix(
  ~ diagnosis + brain_region + sex + age + dataset +
    neu + ast + mic + oli + opc + end,
  data = meta_combined
)
fit_ctrl <- lmFit(expr_combat, design_meta_ctrl)
fit_ctrl <- eBayes(fit_ctrl)
meta_ctrl <- topTable(fit_ctrl, coef = "diagnosisAD",
                      number = Inf, adjust.method = "BH", sort.by = "P")
meta_sig_ctrl <- meta_ctrl[meta_ctrl$adj.P.Val < 0.05 &
                             abs(meta_ctrl$logFC) > 0.5, ]
cat("Meta controlled DEGs:", nrow(meta_sig_ctrl), "\n")

# =============================================================================
# 7. CROSS-COHORT INTERSECTION (Meta-cohort vs GBM)
# =============================================================================

gbm_data <- readRDS(here("data", "processed", "DEG_results.rds"))
gbm_sig_uncorr <- gbm_data$gbm_sig_uncorrected
gbm_sig_ctrl   <- gbm_data$gbm_sig_controlled

common_meta_gbm_uncorr <- intersect(rownames(meta_sig_uncorr),
                                    rownames(gbm_sig_uncorr))
common_meta_gbm_ctrl   <- intersect(rownames(meta_sig_ctrl),
                                    rownames(gbm_sig_ctrl))

cat("\n=== Meta-cohort vs GBM intersection ===\n")
cat("Shared DEGs (uncorrected):", length(common_meta_gbm_uncorr), "\n")
cat("Shared DEGs (controlled):", length(common_meta_gbm_ctrl), "\n")

# Directional analysis
if (length(common_meta_gbm_uncorr) > 0) {
  ad_dirs <- sign(meta_sig_uncorr[common_meta_gbm_uncorr, "logFC"])
  gbm_dirs <- sign(gbm_sig_uncorr[common_meta_gbm_uncorr, "logFC"])
  concordant <- sum(ad_dirs == gbm_dirs)
  cat("Direction concordance:",
      sprintf("%.1f%%", 100 * concordant / length(common_meta_gbm_uncorr)),
      "\n")
  cat("  Both DOWN:", sum(ad_dirs < 0 & gbm_dirs < 0), "\n")
  cat("  Both UP:  ", sum(ad_dirs > 0 & gbm_dirs > 0), "\n")
  cat("  AD DOWN, GBM UP:", sum(ad_dirs < 0 & gbm_dirs > 0), "\n")
  cat("  AD UP, GBM DOWN:", sum(ad_dirs > 0 & gbm_dirs < 0), "\n")
}

# =============================================================================
# 8. SAVE
# =============================================================================

saveRDS(list(expression = expr_combat,
             metadata = meta_combined),
        here("data", "processed", "AD_META_combined.rds"))

saveRDS(list(naive = meta_uncorr,
             controlled = meta_ctrl,
             sig_naive = meta_sig_uncorr,
             sig_controlled = meta_sig_ctrl,
             common_naive_AD_GBM = common_meta_gbm_uncorr,
             common_controlled_AD_GBM = common_meta_gbm_ctrl),
        here("data", "processed", "META_DEG_results.rds"))

write.csv(meta_sig_uncorr,
          here("results", "tables", "AD_META_uncorrected_DEGs.csv"))

cat("\n=== Phase 3.5 (meta-cohort + ComBat + DEG) complete ===\n")
cat("Expected results:\n")
cat("  Meta uncorrected DEGs: 271\n")
cat("  Meta controlled DEGs:  0\n")
cat("  Meta vs GBM shared (uncorrected): 198\n")
cat("  Meta vs GBM shared (controlled):  0\n")