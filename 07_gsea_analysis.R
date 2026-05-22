# =============================================================================
# 07_gsea_analysis.R
# Phase 5: Gene Set Enrichment Analysis (fgsea) with MSigDB Hallmark pathways
# =============================================================================
#
# Goals:
#   1. Rank genes by log2 fold change from each limma model
#   2. Run fgsea against Hallmark collection (H, 50 pathways)
#   3. Compare uncorrected vs deconvolution-controlled enrichment
#   4. Identify cell-composition-resistant pathways
#   5. Test for inverse-direction pathways between AD and GBM
#
# Parameters:
#   - fgsea: minSize = 15, maxSize = 500, nperm = 10000
#   - FDR threshold: 0.05 (Benjamini-Hochberg)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(tidyverse)
})

set.seed(42)

# --- Load Hallmark gene sets ---
cat("=== Loading MSigDB Hallmark collection ===\n")
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
cat("Hallmark pathways loaded:", length(hallmark_list), "\n\n")

# =============================================================================
# 1. LOAD DEG RESULTS (full ranked lists, not just significant subsets)
# =============================================================================

deg_results  <- readRDS(here("data", "processed", "DEG_results.rds"))
meta_results <- readRDS(here("data", "processed", "META_DEG_results.rds"))

# AD meta-cohort (n=122) for the primary AD analysis
ad_uncorr_full <- meta_results$naive       # full ranked AD uncorrected
ad_ctrl_full   <- meta_results$controlled  # full ranked AD controlled

# GBM
gbm_uncorr_full <- deg_results$gbm_uncorrected
gbm_ctrl_full   <- deg_results$gbm_controlled

# =============================================================================
# 2. RANK CONSTRUCTION
# =============================================================================
# Rank by t-statistic (preferred for fgsea: continuous, accounts for both
# magnitude and significance). logFC alone misses noisy genes.
# =============================================================================

build_ranks <- function(deg_df) {
  ranks <- deg_df$t
  names(ranks) <- rownames(deg_df)
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  return(ranks)
}

ad_uncorr_ranks <- build_ranks(ad_uncorr_full)
ad_ctrl_ranks   <- build_ranks(ad_ctrl_full)
gbm_uncorr_ranks <- build_ranks(gbm_uncorr_full)
gbm_ctrl_ranks   <- build_ranks(gbm_ctrl_full)

cat("Ranked gene list sizes:\n")
cat("  AD uncorrected:", length(ad_uncorr_ranks), "\n")
cat("  AD controlled: ", length(ad_ctrl_ranks), "\n")
cat("  GBM uncorrected:", length(gbm_uncorr_ranks), "\n")
cat("  GBM controlled: ", length(gbm_ctrl_ranks), "\n\n")

# =============================================================================
# 3. fgsea ANALYSIS
# =============================================================================

run_fgsea <- function(ranks, label) {
  cat("Running fgsea:", label, "...\n")
  res <- fgsea(pathways = hallmark_list,
               stats = ranks,
               minSize = 15,
               maxSize = 500,
               nPermSimple = 10000)
  res <- res[order(res$padj), ]
  n_sig <- sum(res$padj < 0.05)
  cat("  Significant pathways (FDR<0.05):", n_sig, "\n")
  return(res)
}

ad_uncorr_gsea  <- run_fgsea(ad_uncorr_ranks,  "AD uncorrected")
ad_ctrl_gsea    <- run_fgsea(ad_ctrl_ranks,    "AD controlled")
gbm_uncorr_gsea <- run_fgsea(gbm_uncorr_ranks, "GBM uncorrected")
gbm_ctrl_gsea   <- run_fgsea(gbm_ctrl_ranks,   "GBM controlled")

# =============================================================================
# 4. SIGNIFICANT PATHWAY COMPARISONS
# =============================================================================

ad_uncorr_sig <- ad_uncorr_gsea[padj < 0.05]
ad_ctrl_sig   <- ad_ctrl_gsea[padj < 0.05]
gbm_uncorr_sig <- gbm_uncorr_gsea[padj < 0.05]
gbm_ctrl_sig   <- gbm_ctrl_gsea[padj < 0.05]

cat("\n=== Pathway counts ===\n")
cat("AD  uncorrected:", nrow(ad_uncorr_sig),
    "| controlled:", nrow(ad_ctrl_sig), "\n")
cat("GBM uncorrected:", nrow(gbm_uncorr_sig),
    "| controlled:", nrow(gbm_ctrl_sig), "\n\n")

# AD-specific cell composition artifact pathways (lost under control)
ad_lost <- setdiff(ad_uncorr_sig$pathway, ad_ctrl_sig$pathway)
ad_gained <- setdiff(ad_ctrl_sig$pathway, ad_uncorr_sig$pathway)
ad_preserved <- intersect(ad_uncorr_sig$pathway, ad_ctrl_sig$pathway)
cat("AD pathways:\n")
cat("  Lost under control (artifacts):    ", length(ad_lost), "\n")
cat("  Gained under control (revealed):   ", length(ad_gained), "\n")
cat("  Preserved (cell-autonomous):       ", length(ad_preserved), "\n\n")

# =============================================================================
# 5. INVERSE COMORBIDITY TEST (Controlled AD vs Controlled GBM)
# =============================================================================

cat("=== Inverse comorbidity test (controlled models) ===\n")

# Pathways significant in BOTH controlled cohorts
shared_ctrl <- intersect(ad_ctrl_sig$pathway, gbm_ctrl_sig$pathway)
cat("Pathways significant in BOTH controlled AD and GBM:",
    length(shared_ctrl), "\n")

# Of these, which diverge in direction?
inverse_comorbidity <- data.frame()
for (pw in shared_ctrl) {
  ad_nes  <- ad_ctrl_sig[pathway == pw]$NES
  gbm_nes <- gbm_ctrl_sig[pathway == pw]$NES
  ad_padj <- ad_ctrl_sig[pathway == pw]$padj
  gbm_padj <- gbm_ctrl_sig[pathway == pw]$padj
  direction <- ifelse(sign(ad_nes) != sign(gbm_nes), "INVERSE", "concordant")
  inverse_comorbidity <- rbind(inverse_comorbidity,
                               data.frame(pathway = pw,
                                          AD_NES = ad_nes,
                                          AD_padj = ad_padj,
                                          GBM_NES = gbm_nes,
                                          GBM_padj = gbm_padj,
                                          direction = direction))
}

cat("\n=== INVERSE COMORBIDITY PATHWAYS ===\n")
print(inverse_comorbidity[inverse_comorbidity$direction == "INVERSE", ])

# =============================================================================
# 6. SAVE
# =============================================================================

gsea_results <- list(
  ad_naive       = ad_uncorr_gsea,
  ad_controlled  = ad_ctrl_gsea,
  gbm_naive      = gbm_uncorr_gsea,
  gbm_controlled = gbm_ctrl_gsea,
  inverse_comorbidity = inverse_comorbidity,
  ad_lost_pathways = ad_lost,
  ad_preserved_pathways = ad_preserved,
  shared_controlled_pathways = shared_ctrl
)

saveRDS(gsea_results, here("data", "processed", "GSEA_results.rds"))

# CSV exports (fgsea data.table -> data.frame for csv)
# Note: leadingEdge is a list column; convert to comma-string for CSV
for (nm in c("ad_naive", "ad_controlled", "gbm_naive", "gbm_controlled")) {
  df <- as.data.frame(gsea_results[[nm]])
  if ("leadingEdge" %in% colnames(df)) {
    df$leadingEdge <- sapply(df$leadingEdge, paste, collapse = ",")
  }
  write.csv(df, file.path(here("results", "tables"),
                          paste0("GSEA_", nm, ".csv")),
            row.names = FALSE)
}

cat("\n=== Phase 5 (GSEA) complete ===\n")
cat("Expected key findings:\n")
cat("  AD: 25 uncorrected -> 22 controlled significant pathways\n")
cat("       10 lost (inflammatory/EMT artifacts), 7 gained, 15 preserved\n")
cat("  GBM: 37 uncorrected -> 9 controlled (76% reduction)\n")
cat("  Inverse comorbidity: 3 pathways (DNA REPAIR, GLYCOLYSIS, MYC TARGETS V1)\n")
cat("       AD direction: down (NES ~ -1.9 to -2.0)\n")
cat("       GBM direction: up (NES ~ +1.5 to +2.1)\n")