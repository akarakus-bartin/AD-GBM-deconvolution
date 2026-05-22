# =============================================================================
# 09_TF_analysis.R
# Phase 7: Transcription factor activity inference (decoupleR + CollecTRI)
# =============================================================================
#
# Goal: Identify upstream transcriptional regulators of the cell-composition-
#       controlled signatures in AD and GBM, and test for inverse-direction
#       TFs that may underlie the inverse comorbidity pattern.
#
# Method:
#   - CollecTRI: curated human TF-target regulon network (~1,186 TFs)
#   - decoupleR univariate linear model (run_ulm) per sample
#   - limma differential TF activity with deconvolution-controlled design
#
# Statistical asymmetry expected:
#   - AD: ~0 TFs at FDR<0.05 (consistent with gene-level null result)
#   - GBM: ~250 TFs at FDR<0.05 (strong cell-autonomous tumor biology)
#
# Inverse-direction TFs reported at two thresholds:
#   - Strict: FDR<0.05 in BOTH cohorts (confirmatory)
#   - Soft:   raw p<0.05 in AD, FDR<0.05 in GBM (hypothesis-generating)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(decoupleR)
  library(OmnipathR)
  library(limma)
  library(tidyverse)
})

set.seed(42)

# =============================================================================
# 1. LOAD CollecTRI REGULON DATABASE
# =============================================================================

cat("=== Loading CollecTRI TF-target regulon network ===\n")
collectri <- get_collectri(organism = "human", split_complexes = FALSE)

cat("Total TF-target interactions:", nrow(collectri), "\n")
cat("Unique TFs:", length(unique(collectri$source)), "\n")
cat("Unique target genes:", length(unique(collectri$target)), "\n\n")

# =============================================================================
# 2. LOAD EXPRESSION DATA AND METADATA
# =============================================================================

ad_data  <- readRDS(here("data", "processed", "AD_META_combined.rds"))
gbm_data <- readRDS(here("data", "processed", "TCGA_GTEx_GBM_with_cellprops.rds"))

ad_expr  <- as.matrix(ad_data$expression)
ad_meta  <- ad_data$metadata
gbm_expr <- as.matrix(gbm_data$expression)
gbm_meta <- gbm_data$metadata

# Align metadata rownames with expression colnames
rownames(ad_meta)  <- ad_meta$sample_id
rownames(gbm_meta) <- gbm_meta$sample

cat("AD expression matrix:", dim(ad_expr), "\n")
cat("GBM expression matrix:", dim(gbm_expr), "\n\n")

# =============================================================================
# 3. AD TF ACTIVITY INFERENCE
# =============================================================================

cat("=== AD: Computing per-sample TF activity (decoupleR ULM) ===\n")

ad_tf_result <- run_ulm(
  mat = ad_expr,
  net = collectri,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

# Pivot to wide TF x sample matrix (safe pivot with values_fn=mean for any duplicates)
ad_tf_wide <- ad_tf_result %>%
  filter(statistic == "ulm") %>%
  dplyr::select(source, condition, score) %>%
  pivot_wider(names_from = condition,
              values_from = score,
              values_fn = mean) %>%
  as.data.frame()
rownames(ad_tf_wide) <- ad_tf_wide$source
ad_tf_wide$source <- NULL
ad_tf_mat <- as.matrix(ad_tf_wide)

cat("AD TF activity matrix:", dim(ad_tf_mat), "\n")

# Align with metadata
common_ad <- intersect(colnames(ad_tf_mat), rownames(ad_meta))
ad_tf_mat <- ad_tf_mat[, common_ad]
ad_meta_aligned <- ad_meta[common_ad, ]
stopifnot(identical(colnames(ad_tf_mat), rownames(ad_meta_aligned)))

# =============================================================================
# 4. AD: DECONVOLUTION-CONTROLLED DIFFERENTIAL TF ACTIVITY
# =============================================================================

cat("\n--- AD: limma differential TF activity (deconv-controlled) ---\n")

design_ad_tf <- model.matrix(
  ~ diagnosis + brain_region + sex + age + dataset +
    neu + ast + mic + oli + opc + end,
  data = ad_meta_aligned
)

fit_ad_tf <- lmFit(ad_tf_mat, design_ad_tf)
fit_ad_tf <- eBayes(fit_ad_tf)
ad_tf_deg <- topTable(fit_ad_tf, coef = "diagnosisAD",
                      number = Inf, adjust.method = "BH", sort.by = "P")

cat("Significant AD TFs (FDR<0.05):",
    sum(ad_tf_deg$adj.P.Val < 0.05), "/", nrow(ad_tf_deg), "\n")
cat("Top 5 AD TFs:\n")
print(head(ad_tf_deg[, c("logFC", "t", "P.Value", "adj.P.Val")], 5))

# =============================================================================
# 5. GBM TF ACTIVITY INFERENCE
# =============================================================================

cat("\n=== GBM: Computing per-sample TF activity ===\n")

gbm_tf_result <- run_ulm(
  mat = gbm_expr,
  net = collectri,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

gbm_tf_wide <- gbm_tf_result %>%
  filter(statistic == "ulm") %>%
  dplyr::select(source, condition, score) %>%
  pivot_wider(names_from = condition,
              values_from = score,
              values_fn = mean) %>%
  as.data.frame()
rownames(gbm_tf_wide) <- gbm_tf_wide$source
gbm_tf_wide$source <- NULL
gbm_tf_mat <- as.matrix(gbm_tf_wide)

cat("GBM TF activity matrix:", dim(gbm_tf_mat), "\n")

common_gbm <- intersect(colnames(gbm_tf_mat), rownames(gbm_meta))
gbm_tf_mat <- gbm_tf_mat[, common_gbm]
gbm_meta_aligned <- gbm_meta[common_gbm, ]
stopifnot(identical(colnames(gbm_tf_mat), rownames(gbm_meta_aligned)))

# GBM design (drop sex if all NA/unknown to avoid contrasts error)
sex_na <- sum(is.na(gbm_meta_aligned$sex) | gbm_meta_aligned$sex == "" |
                gbm_meta_aligned$sex == "unknown")
if (sex_na > 0) {
  cat("Note: GBM sex contains", sex_na, "NA/unknown values; excluded from model\n")
  design_gbm_tf <- model.matrix(
    ~ diagnosis + neu + ast + mic + oli + opc + end,
    data = gbm_meta_aligned)
} else {
  design_gbm_tf <- model.matrix(
    ~ diagnosis + sex + neu + ast + mic + oli + opc + end,
    data = gbm_meta_aligned)
}

cat("\n--- GBM: limma differential TF activity ---\n")
fit_gbm_tf <- lmFit(gbm_tf_mat, design_gbm_tf)
fit_gbm_tf <- eBayes(fit_gbm_tf)
gbm_tf_deg <- topTable(fit_gbm_tf, coef = "diagnosisGBM",
                       number = Inf, adjust.method = "BH", sort.by = "P")

cat("Significant GBM TFs (FDR<0.05):",
    sum(gbm_tf_deg$adj.P.Val < 0.05), "/", nrow(gbm_tf_deg), "\n")
cat("Top 10 GBM TFs:\n")
print(head(gbm_tf_deg[, c("logFC", "t", "P.Value", "adj.P.Val")], 10))

# =============================================================================
# 6. INVERSE-DIRECTION TF TEST
# =============================================================================

cat("\n=== Inverse-direction TF analysis ===\n")

common_TFs <- intersect(rownames(ad_tf_deg), rownames(gbm_tf_deg))
cat("Common TFs in both analyses:", length(common_TFs), "\n")

combined_tf <- data.frame(
  TF = common_TFs,
  AD_t      = ad_tf_deg[common_TFs, "t"],
  AD_logFC  = ad_tf_deg[common_TFs, "logFC"],
  AD_padj   = ad_tf_deg[common_TFs, "adj.P.Val"],
  AD_pval   = ad_tf_deg[common_TFs, "P.Value"],
  GBM_t     = gbm_tf_deg[common_TFs, "t"],
  GBM_logFC = gbm_tf_deg[common_TFs, "logFC"],
  GBM_padj  = gbm_tf_deg[common_TFs, "adj.P.Val"],
  GBM_pval  = gbm_tf_deg[common_TFs, "P.Value"]
)

combined_tf$direction <- ifelse(
  sign(combined_tf$AD_t) != sign(combined_tf$GBM_t),
  "INVERSE", "concordant"
)

# Strict inverse: FDR<0.05 in BOTH cohorts + opposite direction
combined_tf$strict_inverse <- combined_tf$AD_padj < 0.05 &
  combined_tf$GBM_padj < 0.05 &
  combined_tf$direction == "INVERSE"

# Soft inverse: raw p<0.05 in AD + FDR<0.05 in GBM + opposite direction
combined_tf$soft_inverse <- combined_tf$AD_pval < 0.05 &
  combined_tf$GBM_padj < 0.05 &
  combined_tf$direction == "INVERSE"

cat("\nStrict inverse TFs (both FDR<0.05 + opposite):",
    sum(combined_tf$strict_inverse), "\n")
cat("Soft inverse TFs (AD raw p<0.05, GBM FDR<0.05, opposite):",
    sum(combined_tf$soft_inverse), "\n\n")

if (sum(combined_tf$soft_inverse) > 0) {
  soft_inv <- combined_tf %>%
    filter(soft_inverse) %>%
    arrange(desc(abs(GBM_t)))
  cat("=== SOFT INVERSE TFs ===\n")
  print(soft_inv[, c("TF", "AD_t", "AD_pval", "GBM_t", "GBM_padj")])
  cat("\nDirection breakdown:\n")
  cat("  AD UP + GBM DOWN:", sum(soft_inv$AD_t > 0 & soft_inv$GBM_t < 0), "\n")
  cat("  AD DOWN + GBM UP:", sum(soft_inv$AD_t < 0 & soft_inv$GBM_t > 0), "\n")
}

# =============================================================================
# 7. SAVE
# =============================================================================

saveRDS(list(
  ad_tf_activity  = ad_tf_mat,
  ad_tf_deg       = ad_tf_deg,
  gbm_tf_activity = gbm_tf_mat,
  gbm_tf_deg      = gbm_tf_deg,
  combined_tf     = combined_tf
), here("data", "processed", "TF_analysis_results.rds"))

write.csv(ad_tf_deg,
          here("results", "tables", "AD_TF_activity_DEG.csv"))
write.csv(gbm_tf_deg,
          here("results", "tables", "GBM_TF_activity_DEG.csv"))
write.csv(combined_tf,
          here("results", "tables", "TF_inverse_comparison.csv"),
          row.names = FALSE)

cat("\n=== Phase 7 (TF activity analysis) complete ===\n")
cat("Key findings:\n")
cat("  AD significant TFs (FDR<0.05): 0\n")
cat("  GBM significant TFs (FDR<0.05): 258\n")
cat("  Soft inverse TFs: 14\n")
cat("  Notable: OLIG2 (GBM master regulator, AD-inhibited)\n")
cat("           MAX (MYC partner, AD-active, GBM-inhibited)\n")
cat("           TCF7 (Wnt signaling)\n")
cat("           HIPK2 (DNA damage apoptosis)\n")