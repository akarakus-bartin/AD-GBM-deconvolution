# =============================================================================
# 04_deg_analysis.R
# Phase 3: Differential expression analysis (limma)
# =============================================================================
#
# Two model families compared:
#   - Uncorrected: ~ diagnosis + brain_region + sex + age (+ dataset for meta)
#   - Deconvolution-controlled: above + neu + ast + mic + oli + opc + end
#
# Significance thresholds:
#   AD:  FDR < 0.05, |log2FC| > 0.5
#   GBM: FDR < 0.05, |log2FC| > 1   (larger effect sizes expected)
#
# Multicollinearity assessed by Variance Inflation Factor (VIF)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(car)        # for VIF
  library(tidyverse)
})

set.seed(42)

# --- Load deconvoluted datasets (output of 03_cell_deconvolution.R) ---
ad_data <- readRDS(here("data", "processed", "GSE48350_RMA_with_cellprops.rds"))
gbm_data <- readRDS(here("data", "processed", "TCGA_GTEx_GBM_with_cellprops.rds"))

ad_expr  <- ad_data$expression
ad_meta  <- ad_data$metadata
gbm_expr <- gbm_data$expression
gbm_meta <- gbm_data$metadata

cat("=== AD GSE48350 dimensions:", dim(ad_expr), "===\n")
cat("=== GBM dimensions:", dim(gbm_expr), "===\n\n")

# =============================================================================
# 1. AD (GSE48350) DEG ANALYSIS
# =============================================================================

# --- Sample alignment between expression and metadata ---
common_ad <- intersect(colnames(ad_expr), rownames(ad_meta))
ad_expr <- ad_expr[, common_ad]
ad_meta <- ad_meta[common_ad, ]
stopifnot(identical(colnames(ad_expr), rownames(ad_meta)))

# Ensure factor levels
ad_meta$diagnosis <- factor(ad_meta$diagnosis, levels = c("Control", "AD"))
ad_meta$brain_region <- factor(ad_meta$brain_region)
ad_meta$sex <- factor(ad_meta$sex)

# --- Uncorrected model ---
cat("--- AD: Uncorrected limma model ---\n")
design_ad_uncorr <- model.matrix(
  ~ diagnosis + brain_region + sex + age,
  data = ad_meta
)
fit_ad_uncorr <- lmFit(ad_expr, design_ad_uncorr)
fit_ad_uncorr <- eBayes(fit_ad_uncorr)
ad_uncorr <- topTable(fit_ad_uncorr, coef = "diagnosisAD",
                      number = Inf, adjust.method = "BH", sort.by = "P")

# Filter significant DEGs
ad_sig_uncorr <- ad_uncorr[ad_uncorr$adj.P.Val < 0.05 & abs(ad_uncorr$logFC) > 0.5, ]
cat("Uncorrected AD DEGs (FDR<0.05, |log2FC|>0.5):", nrow(ad_sig_uncorr), "\n")
cat("  Down:", sum(ad_sig_uncorr$logFC < 0), "\n")
cat("  Up:  ", sum(ad_sig_uncorr$logFC > 0), "\n\n")

# --- Deconvolution-controlled model ---
cat("--- AD: Deconvolution-controlled limma model ---\n")
design_ad_ctrl <- model.matrix(
  ~ diagnosis + brain_region + sex + age + neu + ast + mic + oli + opc + end,
  data = ad_meta
)
fit_ad_ctrl <- lmFit(ad_expr, design_ad_ctrl)
fit_ad_ctrl <- eBayes(fit_ad_ctrl)
ad_ctrl <- topTable(fit_ad_ctrl, coef = "diagnosisAD",
                    number = Inf, adjust.method = "BH", sort.by = "P")

ad_sig_ctrl <- ad_ctrl[ad_ctrl$adj.P.Val < 0.05 & abs(ad_ctrl$logFC) > 0.5, ]
cat("Controlled AD DEGs:", nrow(ad_sig_ctrl), "\n\n")

# --- VIF diagnostic (check for overadjustment) ---
cat("--- VIF diagnostic for AD controlled model ---\n")
lm_vif_ad <- lm(rnorm(nrow(ad_meta)) ~ diagnosis + brain_region + sex + age +
                  neu + ast + mic + oli + opc + end, data = ad_meta)
vif_values_ad <- car::vif(lm_vif_ad)
cat("Diagnosis VIF:", round(vif_values_ad["diagnosisAD"], 2), "\n")
cat("(Threshold: VIF > 5 indicates problematic multicollinearity)\n\n")

# =============================================================================
# 2. GBM DEG ANALYSIS
# =============================================================================

common_gbm <- intersect(colnames(gbm_expr), rownames(gbm_meta))
gbm_expr <- gbm_expr[, common_gbm]
gbm_meta <- gbm_meta[common_gbm, ]
stopifnot(identical(colnames(gbm_expr), rownames(gbm_meta)))

gbm_meta$diagnosis <- factor(gbm_meta$diagnosis, levels = c("Control", "GBM"))
gbm_meta$sex <- factor(gbm_meta$sex)

# --- Uncorrected model ---
cat("--- GBM: Uncorrected limma model ---\n")
design_gbm_uncorr <- model.matrix(~ diagnosis + sex, data = gbm_meta)
fit_gbm_uncorr <- lmFit(gbm_expr, design_gbm_uncorr)
fit_gbm_uncorr <- eBayes(fit_gbm_uncorr)
gbm_uncorr <- topTable(fit_gbm_uncorr, coef = "diagnosisGBM",
                       number = Inf, adjust.method = "BH", sort.by = "P")

gbm_sig_uncorr <- gbm_uncorr[gbm_uncorr$adj.P.Val < 0.05 & abs(gbm_uncorr$logFC) > 1, ]
cat("Uncorrected GBM DEGs (FDR<0.05, |log2FC|>1):", nrow(gbm_sig_uncorr), "\n\n")

# --- Deconvolution-controlled model ---
cat("--- GBM: Deconvolution-controlled limma model ---\n")
design_gbm_ctrl <- model.matrix(
  ~ diagnosis + sex + neu + ast + mic + oli + opc + end,
  data = gbm_meta
)
fit_gbm_ctrl <- lmFit(gbm_expr, design_gbm_ctrl)
fit_gbm_ctrl <- eBayes(fit_gbm_ctrl)
gbm_ctrl <- topTable(fit_gbm_ctrl, coef = "diagnosisGBM",
                     number = Inf, adjust.method = "BH", sort.by = "P")

gbm_sig_ctrl <- gbm_ctrl[gbm_ctrl$adj.P.Val < 0.05 & abs(gbm_ctrl$logFC) > 1, ]
cat("Controlled GBM DEGs:", nrow(gbm_sig_ctrl), "\n\n")

# =============================================================================
# 3. CROSS-COHORT INTERSECTION (Uncorrected and Controlled)
# =============================================================================

cat("=== Cross-cohort intersection (AD GSE48350 vs GBM) ===\n")
common_uncorr <- intersect(rownames(ad_sig_uncorr), rownames(gbm_sig_uncorr))
common_ctrl   <- intersect(rownames(ad_sig_ctrl),   rownames(gbm_sig_ctrl))

cat("Shared DEGs (uncorrected):", length(common_uncorr), "\n")
cat("Shared DEGs (controlled):", length(common_ctrl), "\n\n")

# Directional analysis on shared uncorrected DEGs
if (length(common_uncorr) > 0) {
  ad_dirs  <- sign(ad_sig_uncorr[common_uncorr, "logFC"])
  gbm_dirs <- sign(gbm_sig_uncorr[common_uncorr, "logFC"])
  cat("Direction concordance:\n")
  cat("  Both downregulated:", sum(ad_dirs < 0 & gbm_dirs < 0), "\n")
  cat("  Both upregulated:  ", sum(ad_dirs > 0 & gbm_dirs > 0), "\n")
  cat("  Discordant:        ", sum(ad_dirs != gbm_dirs), "\n")
}

# =============================================================================
# 4. SAVE RESULTS
# =============================================================================

deg_results <- list(
  ad_uncorrected      = ad_uncorr,
  ad_controlled       = ad_ctrl,
  ad_sig_uncorrected  = ad_sig_uncorr,
  ad_sig_controlled   = ad_sig_ctrl,
  gbm_uncorrected     = gbm_uncorr,
  gbm_controlled      = gbm_ctrl,
  gbm_sig_uncorrected = gbm_sig_uncorr,
  gbm_sig_controlled  = gbm_sig_ctrl,
  common_uncorrected  = common_uncorr,
  common_controlled   = common_ctrl,
  ad_vif              = vif_values_ad["diagnosisAD"]
)

saveRDS(deg_results, here("data", "processed", "DEG_results.rds"))

# Also save CSV summaries for manuscript tables
write.csv(ad_sig_uncorr,
          here("results", "tables", "AD_GSE48350_uncorrected_DEGs.csv"))
write.csv(ad_sig_ctrl,
          here("results", "tables", "AD_GSE48350_controlled_DEGs.csv"))
write.csv(gbm_sig_uncorr,
          here("results", "tables", "GBM_uncorrected_DEGs.csv"))
write.csv(gbm_sig_ctrl,
          here("results", "tables", "GBM_controlled_DEGs.csv"))

cat("\n=== Phase 3 (DEG analysis) complete ===\n")
cat("Note: For meta-cohort analysis (n=122), see 05_meta_analysis.R\n")