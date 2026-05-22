# =============================================================================
# 100_supplementary_figures.R
# Supplementary Figures S1-S4 generation (publication-quality, ggplot-based)
# =============================================================================
#
# Figures:
#   S1 - Post-ComBat PCA (dataset effect removed, diagnosis emerges)
#   S2 - AD WGCNA module-trait correlation heatmap
#   S3 - GBM WGCNA module-trait correlation heatmap
#   S4 - AD-GBM module overlap matrix (Fisher's exact test)
#
# All figures use ggplot2 for paper-quality vector output (PDF via cairo).
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(WGCNA)
})

# --- Theme and output directory ---
publication_theme <- theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

fig_sup_dir <- here("results", "figures", "publication")
if (!dir.exists(fig_sup_dir)) dir.create(fig_sup_dir, recursive = TRUE)

# --- Load data ---
ad_meta_data <- readRDS(here("data", "processed", "AD_META_combined.rds"))
wgcna_ad     <- readRDS(here("data", "processed", "WGCNA_AD.rds"))
wgcna_gbm    <- readRDS(here("data", "processed", "WGCNA_GBM.rds"))

cat("=== Generating supplementary figures ===\n\n")

# =============================================================================
# SUPPLEMENTARY FIGURE S1 - Post-ComBat PCA
# =============================================================================

cat("Sup Fig S1: ComBat post-correction PCA...\n")

meta_combined <- ad_meta_data$metadata
expr_combat   <- ad_meta_data$expression

# Use top 2000 variable genes for PCA
combat_var <- apply(expr_combat, 1, var, na.rm = TRUE)
top_var <- order(combat_var, decreasing = TRUE)[1:2000]
pca_post <- prcomp(t(expr_combat[top_var, ]), scale. = TRUE)
ve_post <- pca_post$sdev^2 / sum(pca_post$sdev^2) * 100

pca_df_post <- data.frame(
  PC1 = pca_post$x[, 1],
  PC2 = pca_post$x[, 2],
  dataset   = meta_combined$dataset,
  diagnosis = meta_combined$diagnosis
)

# Compute key correlations for annotation
cor_pc1_dataset <- round(cor(pca_df_post$PC1,
                             as.numeric(factor(pca_df_post$dataset))), 2)
cor_pc1_dx <- round(cor(pca_df_post$PC1,
                        as.numeric(pca_df_post$diagnosis)), 2)

p_sup1 <- ggplot(pca_df_post, aes(x = PC1, y = PC2,
                                  color = dataset, shape = diagnosis)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("GSE48350" = "#1D3557",
                                "GSE36980" = "#E63946")) +
  publication_theme +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = paste0("PC1 x dataset r = ", cor_pc1_dataset,
                          "\nPC1 x diagnosis r = ", cor_pc1_dx),
           size = 3, fontface = "italic") +
  labs(x = paste0("PC1 (", round(ve_post[1], 1), "%)"),
       y = paste0("PC2 (", round(ve_post[2], 1), "%)"),
       color = "Dataset", shape = "Diagnosis")

ggsave(file.path(fig_sup_dir, "FigureS1_ComBat_PCA.pdf"), p_sup1,
       width = 8, height = 6, device = cairo_pdf)
cat("  Saved: FigureS1_ComBat_PCA.pdf\n")

# =============================================================================
# SUPPLEMENTARY FIGURE S2 - AD Module-Trait Correlation Heatmap
# =============================================================================

cat("Sup Fig S2: AD module-trait heatmap...\n")

ad_MEs <- wgcna_ad$MEs
ad_meta_for_wgcna <- wgcna_ad$metadata

ad_traits <- data.frame(
  diagnosis = as.numeric(ad_meta_for_wgcna$diagnosis),
  age       = ad_meta_for_wgcna$age,
  sex       = as.numeric(factor(ad_meta_for_wgcna$sex)),
  hippocampus = as.numeric(ad_meta_for_wgcna$brain_region == "hippocampus"),
  entorhinal  = as.numeric(ad_meta_for_wgcna$brain_region == "entorhinal_cortex"),
  neuron       = ad_meta_for_wgcna$neu,
  astrocyte    = ad_meta_for_wgcna$ast,
  microglia    = ad_meta_for_wgcna$mic,
  oligodendrocyte = ad_meta_for_wgcna$oli,
  OPC          = ad_meta_for_wgcna$opc,
  endothelial  = ad_meta_for_wgcna$end
)

ad_cor  <- cor(ad_MEs, ad_traits, use = "pairwise.complete.obs")
ad_pval <- WGCNA::corPvalueStudent(ad_cor, nrow(ad_MEs))

# Long format for ggplot
ad_long <- expand.grid(Module = rownames(ad_cor), Trait = colnames(ad_cor))
ad_long$r <- as.vector(ad_cor)
ad_long$p <- as.vector(ad_pval)
ad_long$label <- ifelse(
  ad_long$p < 0.001, sprintf("%.2f***", ad_long$r),
  ifelse(ad_long$p < 0.01, sprintf("%.2f**", ad_long$r),
         ifelse(ad_long$p < 0.05, sprintf("%.2f*", ad_long$r),
                sprintf("%.2f", ad_long$r))))

mod_order <- c(paste0("ME", 1:11), "ME0")
mod_order <- mod_order[mod_order %in% unique(ad_long$Module)]
ad_long$Module <- factor(ad_long$Module, levels = rev(mod_order))

trait_order <- c("diagnosis", "age", "sex", "hippocampus", "entorhinal",
                 "neuron", "astrocyte", "microglia", "oligodendrocyte",
                 "OPC", "endothelial")
ad_long$Trait <- factor(ad_long$Trait, levels = trait_order)

p_sup2 <- ggplot(ad_long, aes(x = Trait, y = Module, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  publication_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(x = "Trait", y = "Module")

ggsave(file.path(fig_sup_dir, "FigureS2_AD_module_trait.pdf"), p_sup2,
       width = 10, height = 7, device = cairo_pdf)
cat("  Saved: FigureS2_AD_module_trait.pdf\n")

# =============================================================================
# SUPPLEMENTARY FIGURE S3 - GBM Module-Trait Correlation Heatmap
# =============================================================================

cat("Sup Fig S3: GBM module-trait heatmap...\n")

gbm_MEs <- wgcna_gbm$MEs
gbm_meta_for_wgcna <- wgcna_gbm$metadata
gbm_meta_for_wgcna$sex_num <- ifelse(
  is.na(gbm_meta_for_wgcna$sex) | gbm_meta_for_wgcna$sex == "" |
    gbm_meta_for_wgcna$sex == "unknown",
  NA, as.numeric(factor(gbm_meta_for_wgcna$sex))
)

gbm_traits <- data.frame(
  diagnosis = as.numeric(gbm_meta_for_wgcna$diagnosis),
  sex = gbm_meta_for_wgcna$sex_num,
  tumor       = as.numeric(gbm_meta_for_wgcna$brain_region == "tumor"),
  hippocampus = as.numeric(gbm_meta_for_wgcna$brain_region == "hippocampus"),
  frontal     = as.numeric(gbm_meta_for_wgcna$brain_region == "frontal_cortex"),
  neuron       = gbm_meta_for_wgcna$neu,
  astrocyte    = gbm_meta_for_wgcna$ast,
  microglia    = gbm_meta_for_wgcna$mic,
  oligodendrocyte = gbm_meta_for_wgcna$oli,
  OPC          = gbm_meta_for_wgcna$opc,
  endothelial  = gbm_meta_for_wgcna$end
)

gbm_cor  <- cor(gbm_MEs, gbm_traits, use = "pairwise.complete.obs")
gbm_pval <- WGCNA::corPvalueStudent(gbm_cor, nrow(gbm_MEs))

gbm_long <- expand.grid(Module = rownames(gbm_cor), Trait = colnames(gbm_cor))
gbm_long$r <- as.vector(gbm_cor)
gbm_long$p <- as.vector(gbm_pval)
gbm_long$label <- ifelse(
  gbm_long$p < 0.001, sprintf("%.2f***", gbm_long$r),
  ifelse(gbm_long$p < 0.01, sprintf("%.2f**", gbm_long$r),
         ifelse(gbm_long$p < 0.05, sprintf("%.2f*", gbm_long$r),
                sprintf("%.2f", gbm_long$r))))

gbm_mod_order <- c("ME1", "ME2", "ME3", "ME0")
gbm_long$Module <- factor(gbm_long$Module, levels = rev(gbm_mod_order))

trait_order_gbm <- c("diagnosis", "sex", "tumor", "hippocampus", "frontal",
                     "neuron", "astrocyte", "microglia", "oligodendrocyte",
                     "OPC", "endothelial")
gbm_long$Trait <- factor(gbm_long$Trait, levels = trait_order_gbm)

p_sup3 <- ggplot(gbm_long, aes(x = Trait, y = Module, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  publication_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(x = "Trait", y = "Module")

ggsave(file.path(fig_sup_dir, "FigureS3_GBM_module_trait.pdf"), p_sup3,
       width = 10, height = 4, device = cairo_pdf)
cat("  Saved: FigureS3_GBM_module_trait.pdf\n")

# =============================================================================
# SUPPLEMENTARY FIGURE S4 - AD-GBM Module Overlap Matrix
# =============================================================================

cat("Sup Fig S4: AD-GBM module overlap...\n")

ad_genes_per_module  <- split(colnames(wgcna_ad$expression),
                              wgcna_ad$net$colors)
gbm_genes_per_module <- split(colnames(wgcna_gbm$expression),
                              wgcna_gbm$net$colors)
all_common_genes <- intersect(colnames(wgcna_ad$expression),
                              colnames(wgcna_gbm$expression))
universe_size <- length(all_common_genes)

ad_mods  <- names(ad_genes_per_module)
gbm_mods <- names(gbm_genes_per_module)

overlap_matrix <- matrix(0, nrow = length(ad_mods), ncol = length(gbm_mods),
                         dimnames = list(paste0("AD_M", ad_mods),
                                         paste0("GBM_M", gbm_mods)))
pval_matrix <- overlap_matrix

for (i in seq_along(ad_mods)) {
  ad_mod_genes <- intersect(ad_genes_per_module[[ad_mods[i]]], all_common_genes)
  for (j in seq_along(gbm_mods)) {
    gbm_mod_genes <- intersect(gbm_genes_per_module[[gbm_mods[j]]],
                               all_common_genes)
    a <- length(intersect(ad_mod_genes, gbm_mod_genes))
    b <- length(ad_mod_genes) - a
    c <- length(gbm_mod_genes) - a
    d <- universe_size - a - b - c
    overlap_matrix[i, j] <- a
    if (a > 0 && b >= 0 && c >= 0 && d >= 0) {
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2),
                        alternative = "greater")
      pval_matrix[i, j] <- ft$p.value
    } else {
      pval_matrix[i, j] <- 1
    }
  }
}

overlap_df <- expand.grid(
  AD_module  = rownames(overlap_matrix),
  GBM_module = colnames(overlap_matrix)
)
overlap_df$overlap <- as.vector(overlap_matrix)
overlap_df$log_pval <- as.vector(-log10(pval_matrix + 1e-300))
overlap_df$log_pval[overlap_df$log_pval > 50] <- 50

n_tests <- length(ad_mods) * length(gbm_mods)
bonf_threshold <- 0.05 / n_tests
overlap_df$sig_label <- ifelse(
  as.vector(pval_matrix) < bonf_threshold,
  paste0(overlap_df$overlap, "*"),
  as.character(overlap_df$overlap)
)

# Module ordering by size descending, grey last
ad_sizes <- sapply(ad_genes_per_module, length)
ad_order <- names(sort(ad_sizes, decreasing = TRUE))
ad_order <- c(setdiff(ad_order, "0"), "0")
ad_order_labels <- paste0("AD_M", ad_order)

gbm_sizes <- sapply(gbm_genes_per_module, length)
gbm_order <- names(sort(gbm_sizes, decreasing = TRUE))
gbm_order <- c(setdiff(gbm_order, "0"), "0")
gbm_order_labels <- paste0("GBM_M", gbm_order)

overlap_df$AD_module  <- factor(overlap_df$AD_module,
                                levels = ad_order_labels)
overlap_df$GBM_module <- factor(overlap_df$GBM_module,
                                levels = gbm_order_labels)

p_sup4 <- ggplot(overlap_df, aes(x = GBM_module, y = AD_module, fill = log_pval)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 3.2, fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red",
                       midpoint = 10,
                       name = "-log10(P)\nFisher exact") +
  publication_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(x = "GBM modules", y = "AD modules")

ggsave(file.path(fig_sup_dir, "FigureS4_module_overlap.pdf"), p_sup4,
       width = 9, height = 7, device = cairo_pdf)
cat("  Saved: FigureS4_module_overlap.pdf\n")

cat("\n=== Supplementary figures complete ===\n")
cat("Output directory:", fig_sup_dir, "\n")