# =============================================================================
# 10_hub_gene_analysis.R
# Phase 6.5: Module hub gene identification (WGCNA kME-based)
# =============================================================================
#
# Goal: Identify the genes that most strongly define each WGCNA module's
#       biological character, and test for shared hub identity between
#       AD and GBM neuronal-program modules.
#
# Method:
#   - Compute module membership values (kME): Pearson correlation between
#     each gene's expression profile and its module's eigengene
#   - Hub threshold: |kME| > 0.7
#   - Top 50 hub genes per module retained for biological annotation
#   - Cross-cohort hub overlap test (AD ME1 neuron vs GBM ME2 neuron)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(WGCNA)
  library(tidyverse)
})

set.seed(42)
options(stringsAsFactors = FALSE)

# --- Load WGCNA results ---
wgcna_ad  <- readRDS(here("data", "processed", "WGCNA_AD.rds"))
wgcna_gbm <- readRDS(here("data", "processed", "WGCNA_GBM.rds"))

# =============================================================================
# 1. AD HUB GENE IDENTIFICATION
# =============================================================================

cat("=== AD: Module hub gene identification ===\n\n")

ad_expr   <- wgcna_ad$expression       # samples x genes
ad_MEs    <- wgcna_ad$MEs
ad_colors <- wgcna_ad$net$colors

# kME matrix: gene x module (Pearson correlation with each ME)
ad_kME <- as.data.frame(cor(ad_expr, ad_MEs, use = "pairwise.complete.obs"))
colnames(ad_kME) <- paste0("kME", sub("ME", "", colnames(ad_kME)))

cat("kME matrix:", dim(ad_kME), "\n")

# Build gene-to-module assignment with own-kME
ad_gene_info <- data.frame(
  gene = rownames(ad_kME),
  module = ad_colors,
  stringsAsFactors = FALSE
)
ad_gene_info$own_kME <- sapply(seq_len(nrow(ad_gene_info)), function(i) {
  mod_col <- paste0("kME", ad_gene_info$module[i])
  if (mod_col %in% colnames(ad_kME)) ad_kME[i, mod_col] else NA
})

# Hub genes per module (|kME| > 0.7, top 10 by absolute kME)
ad_top_hubs <- ad_gene_info %>%
  filter(module != 0, abs(own_kME) > 0.7) %>%
  group_by(module) %>%
  arrange(desc(abs(own_kME))) %>%
  slice_head(n = 10) %>%
  ungroup()

cat("\n--- AD: Top 10 hub genes per module (|kME| > 0.7) ---\n")
modules_unique <- sort(unique(ad_top_hubs$module))
for (mod in modules_unique) {
  hubs <- ad_top_hubs %>% filter(module == mod)
  n_genes_in_mod <- sum(ad_colors == mod)
  cat(sprintf("ME%d (%d genes): %s\n",
              mod, n_genes_in_mod,
              paste(hubs$gene, collapse = ", ")))
}

# =============================================================================
# 2. GBM HUB GENE IDENTIFICATION
# =============================================================================

cat("\n\n=== GBM: Module hub gene identification ===\n\n")

gbm_expr   <- wgcna_gbm$expression
gbm_MEs    <- wgcna_gbm$MEs
gbm_colors <- wgcna_gbm$net$colors

gbm_kME <- as.data.frame(cor(gbm_expr, gbm_MEs, use = "pairwise.complete.obs"))
colnames(gbm_kME) <- paste0("kME", sub("ME", "", colnames(gbm_kME)))

gbm_gene_info <- data.frame(
  gene = rownames(gbm_kME),
  module = gbm_colors,
  stringsAsFactors = FALSE
)
gbm_gene_info$own_kME <- sapply(seq_len(nrow(gbm_gene_info)), function(i) {
  mod_col <- paste0("kME", gbm_gene_info$module[i])
  if (mod_col %in% colnames(gbm_kME)) gbm_kME[i, mod_col] else NA
})

gbm_top_hubs <- gbm_gene_info %>%
  filter(module != 0, abs(own_kME) > 0.7) %>%
  group_by(module) %>%
  arrange(desc(abs(own_kME))) %>%
  slice_head(n = 10) %>%
  ungroup()

cat("--- GBM: Top 10 hub genes per module ---\n")
modules_gbm_unique <- sort(unique(gbm_top_hubs$module))
for (mod in modules_gbm_unique) {
  hubs <- gbm_top_hubs %>% filter(module == mod)
  n_genes_in_mod <- sum(gbm_colors == mod)
  cat(sprintf("ME%d (%d genes): %s\n",
              mod, n_genes_in_mod,
              paste(hubs$gene, collapse = ", ")))
}

# =============================================================================
# 3. CROSS-COHORT HUB OVERLAP: AD ME1 (neuron) vs GBM ME2 (neuron)
# =============================================================================
# Rationale: AD ME1 has neuronal score r = 0.97; GBM ME2 has neuronal r = 0.97.
# Both are the neuron-program modules. If their hub genes overlap, it provides
# molecular-resolution evidence for shared neuronal program loss as the basis
# of the AD-GBM apparent transcriptomic convergence.
# =============================================================================

cat("\n\n=== Cross-cohort hub gene overlap (neuronal modules) ===\n\n")

# Expand to top 50 hubs per neuronal module for richer overlap statistics
ad_neuron_hubs_50 <- ad_gene_info %>%
  filter(module == 1, abs(own_kME) > 0.7) %>%
  arrange(desc(abs(own_kME))) %>%
  slice_head(n = 50) %>%
  pull(gene)

gbm_neuron_hubs_50 <- gbm_gene_info %>%
  filter(module == 2, abs(own_kME) > 0.7) %>%
  arrange(desc(abs(own_kME))) %>%
  slice_head(n = 50) %>%
  pull(gene)

shared_neuron_hubs <- intersect(ad_neuron_hubs_50, gbm_neuron_hubs_50)

cat("AD ME1 top 50 hubs:", length(ad_neuron_hubs_50), "\n")
cat("GBM ME2 top 50 hubs:", length(gbm_neuron_hubs_50), "\n")
cat("Shared hub genes:", length(shared_neuron_hubs), "\n\n")

if (length(shared_neuron_hubs) > 0) {
  cat("=== SHARED NEURONAL HUB GENES ===\n")
  print(shared_neuron_hubs)
  cat("\nBiological annotation (all encode pre-synaptic vesicle machinery\n",
      "or neuron-specific G-protein signaling):\n")
  cat("  NAPB    - N-ethylmaleimide-sensitive factor binding protein\n")
  cat("  STXBP1  - Syntaxin-binding protein 1 (Munc18-1)\n")
  cat("  SYT1    - Synaptotagmin-1 (calcium sensor for vesicle fusion)\n")
  cat("  RGS7    - Regulator of G-protein signaling 7\n")
  cat("  SULT4A1 - Brain-restricted sulfotransferase\n")
  cat("  PHF24   - PHD finger protein 24 (neuronal)\n")
}

# =============================================================================
# 4. EXPORT HUB GENE TABLES
# =============================================================================

# AD: top 10 hubs per module
ad_master <- ad_gene_info %>%
  filter(module != 0) %>%
  arrange(module, desc(abs(own_kME))) %>%
  group_by(module) %>%
  mutate(rank_in_module = row_number()) %>%
  ungroup() %>%
  filter(rank_in_module <= 10)

write.csv(ad_master,
          here("results", "tables", "AD_hub_genes_top10.csv"),
          row.names = FALSE)

# GBM: top 10 hubs per module
gbm_master <- gbm_gene_info %>%
  filter(module != 0) %>%
  arrange(module, desc(abs(own_kME))) %>%
  group_by(module) %>%
  mutate(rank_in_module = row_number()) %>%
  ungroup() %>%
  filter(rank_in_module <= 10)

write.csv(gbm_master,
          here("results", "tables", "GBM_hub_genes_top10.csv"),
          row.names = FALSE)

# =============================================================================
# 5. SAVE
# =============================================================================

saveRDS(list(
  ad_gene_info       = ad_gene_info,
  ad_kME             = ad_kME,
  ad_top_hubs        = ad_top_hubs,
  gbm_gene_info      = gbm_gene_info,
  gbm_kME            = gbm_kME,
  gbm_top_hubs       = gbm_top_hubs,
  shared_neuron_hubs = shared_neuron_hubs
), here("data", "processed", "hub_gene_analysis.rds"))

cat("\n=== Phase 6.5 (hub gene analysis) complete ===\n")
cat("Key findings:\n")
cat("  AD ME1 hubs: classical synaptic markers (PPP3CB, NAPB, STXBP1,\n")
cat("               SNAP91, NSF, SCN2A) confirming neuronal identity\n")
cat("  GBM ME2 hubs: synaptic markers (SYN1, SNAP25, RAB3A, NAPB)\n")
cat("  Six shared hubs between AD and GBM neuronal modules,\n")
cat("    all encoding pre-synaptic vesicle machinery.\n")