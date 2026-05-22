# =============================================================================
# Publication-Ready Figures for AD-GBM Comparative Transcriptomics
# All in English, paper-quality, cairo_pdf for cross-platform compatibility
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(WGCNA)
  library(pheatmap)
  library(RColorBrewer)
  library(fgsea)
})

# Common theme for all figures
publication_theme <- theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Output directory
fig_pub_dir <- here("results", "figures", "publication")
if (!dir.exists(fig_pub_dir)) dir.create(fig_pub_dir, recursive = TRUE)

cat("=== Publication figures being generated ===\n\n")

# =============================================================================
# Load all data
# =============================================================================
ad_data  <- readRDS(here("data", "processed", "AD_META_combined.rds"))
gbm_data <- readRDS(here("data", "processed", "TCGA_GTEx_GBM_with_cellprops.rds"))
deg_results <- readRDS(here("data", "processed", "DEG_results.rds"))
meta_deg_results <- readRDS(here("data", "processed", "META_DEG_results.rds"))
wgcna_ad <- readRDS(here("data", "processed", "WGCNA_AD.rds"))
wgcna_gbm <- readRDS(here("data", "processed", "WGCNA_GBM.rds"))
gsea_results <- readRDS(here("data", "processed", "GSEA_results.rds"))

# =============================================================================
# FIGURE 1: Cell-Type Composition Differences (Key Finding)
# =============================================================================

cat("Figure 1: Cell-type composition...\n")

# Combined data
ad_meta <- ad_data$metadata
ad_long <- ad_meta %>%
  dplyr::select(diagnosis, neu, ast, mic, oli, opc, end) %>%
  pivot_longer(cols = -diagnosis, names_to = "cell_type", values_to = "score") %>%
  mutate(cohort = "AD (Meta, n=122)")

gbm_meta <- gbm_data$metadata
gbm_long <- gbm_meta %>%
  dplyr::select(diagnosis, neu, ast, mic, oli, opc, end) %>%
  pivot_longer(cols = -diagnosis, names_to = "cell_type", values_to = "score") %>%
  mutate(cohort = "GBM (TCGA + GTEx, n=339)")

combined_cells <- bind_rows(ad_long, gbm_long)
combined_cells$cell_type <- factor(
  combined_cells$cell_type,
  levels = c("neu", "ast", "mic", "oli", "opc", "end"),
  labels = c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Endothelial")
)
combined_cells$diagnosis_label <- ifelse(combined_cells$diagnosis == "Control", "Control", "Disease")
combined_cells$diagnosis_label <- factor(combined_cells$diagnosis_label, levels = c("Control", "Disease"))

p_fig1 <- ggplot(combined_cells, aes(x = diagnosis_label, y = score, fill = diagnosis_label)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.6) +
  geom_jitter(width = 0.18, size = 0.3, alpha = 0.3) +
  facet_grid(cohort ~ cell_type, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#457B9D", "Disease" = "#E63946")) +
  publication_theme +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 9)
  ) +
  labs(
    title = "Cell-type composition differences are dramatic in both AD and GBM",
    x = "", y = "BRETIGEA cell-type score"
  )

ggsave(file.path(fig_pub_dir, "Figure1_cell_composition.pdf"), p_fig1,
       width = 11, height = 6, device = cairo_pdf)
cat("  Saved: Figure1_cell_composition.pdf\n")

# =============================================================================
# FIGURE 2: Neuron Loss as Central Finding
# =============================================================================

cat("Figure 2: Neuron loss key finding...\n")

neuron_only <- combined_cells %>% filter(cell_type == "Neuron")

p_fig2 <- ggplot(neuron_only, aes(x = diagnosis_label, y = score, fill = diagnosis_label)) +
  geom_violin(alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.2, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.4) +
  facet_wrap(~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#457B9D", "Disease" = "#E63946")) +
  publication_theme +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 11)) +
  labs(
    title = "Dramatic neuron loss in both AD and GBM tissues",
    subtitle = "Underlies apparent 'shared synaptic signature' in naive analyses",
    x = "", y = "Neuronal cell-type score"
  )

ggsave(file.path(fig_pub_dir, "Figure2_neuron_loss.pdf"), p_fig2,
       width = 8, height = 5, device = cairo_pdf)
cat("  Saved: Figure2_neuron_loss.pdf\n")

# =============================================================================
# FIGURE 3: DEG counts — naive vs deconvolution-controlled
# =============================================================================

cat("Figure 3: DEG comparison...\n")

deg_summary <- data.frame(
  Cohort = c("AD (GSE48350, n=75)", "AD (GSE48350, n=75)",
             "AD (Meta, n=122)", "AD (Meta, n=122)",
             "GBM (TCGA+GTEx, n=339)", "GBM (TCGA+GTEx, n=339)"),
  Model = c("Naive", "Controlled", "Naive", "Controlled", "Naive", "Controlled"),
  DEG_count = c(489, 0, 271, 0, 11082, 3788)
)
deg_summary$Model <- factor(deg_summary$Model, levels = c("Naive", "Controlled"))
deg_summary$Cohort <- factor(deg_summary$Cohort, levels = unique(deg_summary$Cohort))

p_fig3 <- ggplot(deg_summary, aes(x = Cohort, y = DEG_count + 1, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = DEG_count),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  publication_theme +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  labs(
    title = "Differentially expressed genes: naive vs deconvolution-controlled",
    subtitle = "AD signature collapses entirely; GBM retains cell-autonomous tumor biology",
    y = "Significant DEGs (FDR<0.05, log10 scale)", x = ""
  )

ggsave(file.path(fig_pub_dir, "Figure3_DEG_comparison.pdf"), p_fig3,
       width = 9, height = 6, device = cairo_pdf)
cat("  Saved: Figure3_DEG_comparison.pdf\n")

# =============================================================================
# FIGURE 4: Central finding — AD-GBM intersection
# =============================================================================

cat("Figure 4: AD-GBM intersection (central finding)...\n")

central_data <- data.frame(
  Question = c("Naive AD", "Naive GBM", "Shared (Naive)",
               "Controlled AD", "Controlled GBM", "Shared (Controlled)"),
  Count = c(271, 11082, 322, 0, 3788, 0),
  Model = c("Naive", "Naive", "Naive", "Controlled", "Controlled", "Controlled")
)
central_data$Question <- factor(central_data$Question, levels = central_data$Question)
central_data$Model <- factor(central_data$Model, levels = c("Naive", "Controlled"))

p_fig4 <- ggplot(central_data, aes(x = Question, y = Count + 1, fill = Model)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.4, size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  publication_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(
    title = "Central finding: apparent AD-GBM convergence disappears under deconvolution",
    subtitle = "322 'shared' DEGs in naive models reduce to 0 after cell composition control",
    y = "Significant DEGs (log10 scale + 1)", x = ""
  )

ggsave(file.path(fig_pub_dir, "Figure4_central_finding.pdf"), p_fig4,
       width = 10, height = 6, device = cairo_pdf)
cat("  Saved: Figure4_central_finding.pdf\n")

# =============================================================================
# FIGURE 5: Inverse Comorbidity (KEY)
# =============================================================================

cat("Figure 5: Inverse comorbidity (key finding)...\n")

inverse_pathways <- c("HALLMARK_DNA_REPAIR", "HALLMARK_GLYCOLYSIS", "HALLMARK_MYC_TARGETS_V1")

inverse_df <- rbind(
  gsea_results$ad_controlled[pathway %in% inverse_pathways, 
                             .(pathway, NES, padj, Cohort = "AD")],
  gsea_results$gbm_controlled[pathway %in% inverse_pathways, 
                              .(pathway, NES, padj, Cohort = "GBM")]
)
inverse_df$pathway_clean <- gsub("HALLMARK_", "", inverse_df$pathway)
inverse_df$pathway_clean <- gsub("_", " ", inverse_df$pathway_clean)

p_fig5 <- ggplot(inverse_df, aes(x = pathway_clean, y = NES, fill = Cohort)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("p=%.0e", padj)),
            position = position_dodge(width = 0.7),
            vjust = ifelse(inverse_df$NES > 0, -0.5, 1.5),
            size = 3.3) +
  scale_fill_manual(values = c("AD" = "#1D3557", "GBM" = "#E63946")) +
  publication_theme +
  theme(axis.text.x = element_text(face = "bold", size = 10)) +
  labs(
    title = "Cell-composition-resistant pathways diverge in opposite directions",
    subtitle = "Molecular evidence for AD-cancer inverse comorbidity",
    x = "Hallmark pathway", y = "Normalized Enrichment Score (NES)"
  )

ggsave(file.path(fig_pub_dir, "Figure5_inverse_comorbidity.pdf"), p_fig5,
       width = 9, height = 5.5, device = cairo_pdf)
cat("  Saved: Figure5_inverse_comorbidity.pdf\n")

# =============================================================================
# FIGURE 6: AD GSEA — naive vs controlled (all pathways)
# =============================================================================

cat("Figure 6: AD GSEA naive vs controlled...\n")

ad_pathways_all <- unique(c(
  gsea_results$ad_naive[padj < 0.05]$pathway,
  gsea_results$ad_controlled[padj < 0.05]$pathway
))

ad_gsea_df <- rbind(
  gsea_results$ad_naive[pathway %in% ad_pathways_all, 
                        .(pathway, NES, padj, Model = "Naive")],
  gsea_results$ad_controlled[pathway %in% ad_pathways_all,
                             .(pathway, NES, padj, Model = "Controlled")]
)
ad_gsea_df$Model <- factor(ad_gsea_df$Model, levels = c("Naive", "Controlled"))
ad_gsea_df$pathway_clean <- gsub("HALLMARK_", "", ad_gsea_df$pathway)
ad_gsea_df$pathway_clean <- gsub("_", " ", ad_gsea_df$pathway_clean)
ad_gsea_df$significant <- ad_gsea_df$padj < 0.05

# Sort by max abs NES across models
order_ad <- ad_gsea_df %>%
  group_by(pathway_clean) %>%
  summarise(max_abs = max(abs(NES), na.rm = TRUE)) %>%
  arrange(max_abs) %>%
  pull(pathway_clean)
ad_gsea_df$pathway_clean <- factor(ad_gsea_df$pathway_clean, levels = order_ad)

p_fig6 <- ggplot(ad_gsea_df, aes(x = NES, y = pathway_clean, fill = Model, alpha = significant)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  publication_theme +
  labs(
    title = "AD Hallmark pathway enrichment: cell composition versus cell-autonomous signal",
    subtitle = "Inflammatory and EMT pathways disappear under control; metabolic dysfunction persists",
    x = "Normalized Enrichment Score (NES)", y = ""
  )

ggsave(file.path(fig_pub_dir, "Figure6_AD_GSEA.pdf"), p_fig6,
       width = 11, height = 8, device = cairo_pdf)
cat("  Saved: Figure6_AD_GSEA.pdf\n")

# =============================================================================
# FIGURE 7: GBM GSEA — naive vs controlled
# =============================================================================

cat("Figure 7: GBM GSEA naive vs controlled...\n")

gbm_pathways_all <- unique(c(
  gsea_results$gbm_naive[padj < 0.05][order(padj)][1:20]$pathway,
  gsea_results$gbm_controlled[padj < 0.05]$pathway
))

gbm_gsea_df <- rbind(
  gsea_results$gbm_naive[pathway %in% gbm_pathways_all,
                         .(pathway, NES, padj, Model = "Naive")],
  gsea_results$gbm_controlled[pathway %in% gbm_pathways_all,
                              .(pathway, NES, padj, Model = "Controlled")]
)
gbm_gsea_df$Model <- factor(gbm_gsea_df$Model, levels = c("Naive", "Controlled"))
gbm_gsea_df$pathway_clean <- gsub("HALLMARK_", "", gbm_gsea_df$pathway)
gbm_gsea_df$pathway_clean <- gsub("_", " ", gbm_gsea_df$pathway_clean)
gbm_gsea_df$significant <- gbm_gsea_df$padj < 0.05

order_gbm <- gbm_gsea_df %>%
  group_by(pathway_clean) %>%
  summarise(max_abs = max(abs(NES), na.rm = TRUE)) %>%
  arrange(max_abs) %>%
  pull(pathway_clean)
gbm_gsea_df$pathway_clean <- factor(gbm_gsea_df$pathway_clean, levels = order_gbm)

p_fig7 <- ggplot(gbm_gsea_df, aes(x = NES, y = pathway_clean, fill = Model, alpha = significant)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  publication_theme +
  labs(
    title = "GBM Hallmark pathway enrichment: 76% reduction after cell composition control",
    subtitle = "Surviving pathways reflect bona fide tumor cell proliferation programs",
    x = "Normalized Enrichment Score (NES)", y = ""
  )

ggsave(file.path(fig_pub_dir, "Figure7_GBM_GSEA.pdf"), p_fig7,
       width = 11, height = 8, device = cairo_pdf)
cat("  Saved: Figure7_GBM_GSEA.pdf\n")

cat("\n=== All publication figures generated ===\n")
cat("Output directory:", fig_pub_dir, "\n")



# =============================================================================
# Publication figures — title/subtitle KALDIRILMIS versiyon
# Dergi submission standardına uygun (caption ayri verilir)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggplot2)
  library(fgsea)
})

publication_theme <- theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

fig_pub_dir <- here("results", "figures", "publication")
if (!dir.exists(fig_pub_dir)) dir.create(fig_pub_dir, recursive = TRUE)

# Load data
ad_data  <- readRDS(here("data", "processed", "AD_META_combined.rds"))
gbm_data <- readRDS(here("data", "processed", "TCGA_GTEx_GBM_with_cellprops.rds"))
gsea_results <- readRDS(here("data", "processed", "GSEA_results.rds"))

cat("=== Regenerating figures without titles ===\n\n")

# --- FIGURE 1: Cell composition ---
ad_long <- ad_data$metadata %>%
  dplyr::select(diagnosis, neu, ast, mic, oli, opc, end) %>%
  pivot_longer(cols = -diagnosis, names_to = "cell_type", values_to = "score") %>%
  mutate(cohort = "AD (Meta, n=122)")

gbm_long <- gbm_data$metadata %>%
  dplyr::select(diagnosis, neu, ast, mic, oli, opc, end) %>%
  pivot_longer(cols = -diagnosis, names_to = "cell_type", values_to = "score") %>%
  mutate(cohort = "GBM (TCGA + GTEx, n=339)")

combined_cells <- bind_rows(ad_long, gbm_long)
combined_cells$cell_type <- factor(
  combined_cells$cell_type,
  levels = c("neu", "ast", "mic", "oli", "opc", "end"),
  labels = c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Endothelial")
)
combined_cells$diagnosis_label <- factor(
  ifelse(combined_cells$diagnosis == "Control", "Control", "Disease"),
  levels = c("Control", "Disease")
)

p_fig1 <- ggplot(combined_cells, aes(x = diagnosis_label, y = score, fill = diagnosis_label)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.6) +
  geom_jitter(width = 0.18, size = 0.3, alpha = 0.3) +
  facet_grid(cohort ~ cell_type, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#457B9D", "Disease" = "#E63946")) +
  publication_theme +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 9)) +
  labs(x = "", y = "BRETIGEA cell-type score")

ggsave(file.path(fig_pub_dir, "Figure1_cell_composition.pdf"), p_fig1,
       width = 11, height = 6, device = cairo_pdf)
cat("Figure1 done\n")

# --- FIGURE 2: Neuron loss ---
neuron_only <- combined_cells %>% filter(cell_type == "Neuron")

p_fig2 <- ggplot(neuron_only, aes(x = diagnosis_label, y = score, fill = diagnosis_label)) +
  geom_violin(alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.2, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.4) +
  facet_wrap(~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#457B9D", "Disease" = "#E63946")) +
  publication_theme +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 11)) +
  labs(x = "", y = "Neuronal cell-type score")

ggsave(file.path(fig_pub_dir, "Figure2_neuron_loss.pdf"), p_fig2,
       width = 8, height = 5, device = cairo_pdf)
cat("Figure2 done\n")

# --- FIGURE 3: DEG comparison ---
deg_summary <- data.frame(
  Cohort = factor(c("AD (GSE48350, n=75)", "AD (GSE48350, n=75)",
                    "AD (Meta, n=122)", "AD (Meta, n=122)",
                    "GBM (TCGA+GTEx, n=339)", "GBM (TCGA+GTEx, n=339)"),
                  levels = c("AD (GSE48350, n=75)", "AD (Meta, n=122)", "GBM (TCGA+GTEx, n=339)")),
  Model = factor(c("Naive", "Controlled", "Naive", "Controlled", "Naive", "Controlled"),
                 levels = c("Naive", "Controlled")),
  DEG_count = c(489, 0, 271, 0, 11082, 3788)
)

p_fig3 <- ggplot(deg_summary, aes(x = Cohort, y = DEG_count + 1, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = DEG_count),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  publication_theme +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(y = "Significant DEGs (FDR<0.05, log10 scale)", x = "")

ggsave(file.path(fig_pub_dir, "Figure3_DEG_comparison.pdf"), p_fig3,
       width = 9, height = 6, device = cairo_pdf)
cat("Figure3 done\n")


# --- FIGURE 4: Central finding (UPDATED: meta-cohort, n=122 vs n=339) ---
# Manuscript v2 ile tutarli — 198 shared DEGs (AD meta vs GBM)
central_data <- data.frame(
  Question = factor(c("Uncorrected AD\n(Meta, n=122)", 
                      "Uncorrected GBM\n(n=339)", 
                      "Shared\n(Uncorrected)",
                      "Controlled AD\n(Meta, n=122)", 
                      "Controlled GBM\n(n=339)", 
                      "Shared\n(Controlled)"),
                    levels = c("Uncorrected AD\n(Meta, n=122)", 
                               "Uncorrected GBM\n(n=339)", 
                               "Shared\n(Uncorrected)",
                               "Controlled AD\n(Meta, n=122)", 
                               "Controlled GBM\n(n=339)", 
                               "Shared\n(Controlled)")),
  Count = c(271, 11082, 198, 0, 3788, 0),  # 322 -> 198
  Model = factor(c("Uncorrected", "Uncorrected", "Uncorrected", 
                   "Controlled", "Controlled", "Controlled"),
                 levels = c("Uncorrected", "Controlled"))
)

p_fig4 <- ggplot(central_data, aes(x = Question, y = Count + 1, fill = Model)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.4, size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c("Uncorrected" = "#E63946", "Controlled" = "#457B9D")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  publication_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
        axis.title.y = element_text(size = 10, margin = margin(r = 8))) +
  labs(y = "Number of significant DEGs", x = "")

ggsave(file.path(fig_pub_dir, "Figure4_central_finding.pdf"), p_fig4,
       width = 10, height = 6, device = cairo_pdf)
cat("Figure4 updated (198 shared, meta-cohort)\n")

# --- FIGURE 5: Inverse comorbidity ---
inverse_pathways <- c("HALLMARK_DNA_REPAIR", "HALLMARK_GLYCOLYSIS", "HALLMARK_MYC_TARGETS_V1")

inverse_df <- rbind(
  gsea_results$ad_controlled[pathway %in% inverse_pathways, 
                             .(pathway, NES, padj, Cohort = "AD")],
  gsea_results$gbm_controlled[pathway %in% inverse_pathways, 
                              .(pathway, NES, padj, Cohort = "GBM")]
)
inverse_df$pathway_clean <- gsub("HALLMARK_", "", inverse_df$pathway)
inverse_df$pathway_clean <- gsub("_", " ", inverse_df$pathway_clean)

p_fig5 <- ggplot(inverse_df, aes(x = pathway_clean, y = NES, fill = Cohort)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("p=%.0e", padj)),
            position = position_dodge(width = 0.7),
            vjust = ifelse(inverse_df$NES > 0, -0.5, 1.5),
            size = 3.3) +
  scale_fill_manual(values = c("AD" = "#1D3557", "GBM" = "#E63946")) +
  publication_theme +
  theme(axis.text.x = element_text(face = "bold", size = 10)) +
  labs(x = "Hallmark pathway", y = "Normalized Enrichment Score (NES)")

ggsave(file.path(fig_pub_dir, "Figure5_inverse_comorbidity.pdf"), p_fig5,
       width = 9, height = 5.5, device = cairo_pdf)
cat("Figure5 done\n")

# --- FIGURE 6: AD GSEA ---
ad_pathways_all <- unique(c(
  gsea_results$ad_naive[padj < 0.05]$pathway,
  gsea_results$ad_controlled[padj < 0.05]$pathway
))

ad_gsea_df <- rbind(
  gsea_results$ad_naive[pathway %in% ad_pathways_all, 
                        .(pathway, NES, padj, Model = "Naive")],
  gsea_results$ad_controlled[pathway %in% ad_pathways_all,
                             .(pathway, NES, padj, Model = "Controlled")]
)
ad_gsea_df$Model <- factor(ad_gsea_df$Model, levels = c("Naive", "Controlled"))
ad_gsea_df$pathway_clean <- gsub("HALLMARK_", "", ad_gsea_df$pathway)
ad_gsea_df$pathway_clean <- gsub("_", " ", ad_gsea_df$pathway_clean)
ad_gsea_df$significant <- ad_gsea_df$padj < 0.05

order_ad <- ad_gsea_df %>%
  group_by(pathway_clean) %>%
  summarise(max_abs = max(abs(NES), na.rm = TRUE)) %>%
  arrange(max_abs) %>%
  pull(pathway_clean)
ad_gsea_df$pathway_clean <- factor(ad_gsea_df$pathway_clean, levels = order_ad)

p_fig6 <- ggplot(ad_gsea_df, aes(x = NES, y = pathway_clean, fill = Model, alpha = significant)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  publication_theme +
  labs(x = "Normalized Enrichment Score (NES)", y = "")

ggsave(file.path(fig_pub_dir, "Figure6_AD_GSEA.pdf"), p_fig6,
       width = 11, height = 8, device = cairo_pdf)
cat("Figure6 done\n")

# --- FIGURE 7: GBM GSEA ---
gbm_pathways_all <- unique(c(
  gsea_results$gbm_naive[padj < 0.05][order(padj)][1:20]$pathway,
  gsea_results$gbm_controlled[padj < 0.05]$pathway
))

gbm_gsea_df <- rbind(
  gsea_results$gbm_naive[pathway %in% gbm_pathways_all,
                         .(pathway, NES, padj, Model = "Naive")],
  gsea_results$gbm_controlled[pathway %in% gbm_pathways_all,
                              .(pathway, NES, padj, Model = "Controlled")]
)
gbm_gsea_df$Model <- factor(gbm_gsea_df$Model, levels = c("Naive", "Controlled"))
gbm_gsea_df$pathway_clean <- gsub("HALLMARK_", "", gbm_gsea_df$pathway)
gbm_gsea_df$pathway_clean <- gsub("_", " ", gbm_gsea_df$pathway_clean)
gbm_gsea_df$significant <- gbm_gsea_df$padj < 0.05

order_gbm <- gbm_gsea_df %>%
  group_by(pathway_clean) %>%
  summarise(max_abs = max(abs(NES), na.rm = TRUE)) %>%
  arrange(max_abs) %>%
  pull(pathway_clean)
gbm_gsea_df$pathway_clean <- factor(gbm_gsea_df$pathway_clean, levels = order_gbm)

p_fig7 <- ggplot(gbm_gsea_df, aes(x = NES, y = pathway_clean, fill = Model, alpha = significant)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  publication_theme +
  labs(x = "Normalized Enrichment Score (NES)", y = "")

ggsave(file.path(fig_pub_dir, "Figure7_GBM_GSEA.pdf"), p_fig7,
       width = 11, height = 8, device = cairo_pdf)
cat("Figure7 done\n")

cat("\n=== All figures regenerated without titles ===\n")












# =============================================================================
# Figure 1 ve Figure 4 — minor fixes
# =============================================================================

# --- FIGURE 1 fix: y-ekseni daha temiz, facet scales ayri ---
p_fig1 <- ggplot(combined_cells, aes(x = diagnosis_label, y = score, fill = diagnosis_label)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.6) +
  geom_jitter(width = 0.18, size = 0.3, alpha = 0.3) +
  facet_grid(cohort ~ cell_type, scales = "free", switch = "y") +  # free hem x hem y
  scale_fill_manual(values = c("Control" = "#457B9D", "Disease" = "#E63946")) +
  scale_y_continuous(n.breaks = 4) +  # daha az break, daha temiz
  publication_theme +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 9),
        strip.placement = "outside",
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 8)) +
  labs(x = "", y = "BRETIGEA cell-type score")

ggsave(file.path(fig_pub_dir, "Figure1_cell_composition.pdf"), p_fig1,
       width = 12, height = 6, device = cairo_pdf)
cat("Figure1 fixed\n")

# --- FIGURE 4 fix: y-axis text tek satir ---
p_fig4 <- ggplot(central_data, aes(x = Question, y = Count + 1, fill = Model)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.4, size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  publication_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(y = "Significant DEGs (log10 + 1)", x = "")

ggsave(file.path(fig_pub_dir, "Figure4_central_finding.pdf"), p_fig4,
       width = 10, height = 6, device = cairo_pdf)
cat("Figure4 fixed\n")

# --- FIGURE 3 fix: ayni problem var ---
p_fig3 <- ggplot(deg_summary, aes(x = Cohort, y = DEG_count + 1, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = DEG_count),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  publication_theme +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(y = "Significant DEGs (FDR < 0.05)", x = "")

ggsave(file.path(fig_pub_dir, "Figure3_DEG_comparison.pdf"), p_fig3,
       width = 9, height = 6, device = cairo_pdf)
cat("Figure3 fixed\n")







# =============================================================================
# Figure 1 ve Figure 3 — final fixes
# =============================================================================

# --- FIGURE 1: facet_wrap kullan, daha temiz layout ---
# Sorun: facet_grid ile iki cohort yan yana × 6 cell type, 12 panel,
# y-ekseninde sayilar siktigi icin ust uste biniyor.
# Cozum: Her cohort ayri panel (faceted by cohort and cell), 
# y-axis breaks acik sekilde belirtildi

p_fig1 <- ggplot(combined_cells, aes(x = diagnosis_label, y = score, fill = diagnosis_label)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.6) +
  geom_jitter(width = 0.18, size = 0.3, alpha = 0.3) +
  facet_grid(cohort ~ cell_type, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#457B9D", "Disease" = "#E63946")) +
  publication_theme +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        panel.spacing.y = unit(1, "lines"),  # Cohort'lar arasında dikey boşluk
        panel.spacing.x = unit(0.5, "lines")) +
  labs(x = "", y = "BRETIGEA cell-type score")

ggsave(file.path(fig_pub_dir, "Figure1_cell_composition.pdf"), p_fig1,
       width = 13, height = 7, device = cairo_pdf)
cat("Figure1 fixed (wider canvas, more panel spacing)\n")

# --- FIGURE 3: y-ekseni baslik kisa, etiket ayri satira gecmesin ---
# Sorun: "Significant DEGs (FDR < 0.05)" cogu cizicide iki satira bolunuyor
# Cozum: Etiketi kisalt ve y-axis title size'ini ayarla

p_fig3 <- ggplot(deg_summary, aes(x = Cohort, y = DEG_count + 1, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = DEG_count),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  publication_theme +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1),
        axis.title.y = element_text(size = 10, margin = margin(r = 8))) +
  labs(y = "Number of significant DEGs", x = "")

ggsave(file.path(fig_pub_dir, "Figure3_DEG_comparison.pdf"), p_fig3,
       width = 9, height = 6, device = cairo_pdf)
cat("Figure3 fixed (cleaner y-axis label)\n")

# --- FIGURE 4: Ayni y-axis problemi olabilir, kontrol icin yeniden ---
p_fig4 <- ggplot(central_data, aes(x = Question, y = Count + 1, fill = Model)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.4, size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c("Naive" = "#E63946", "Controlled" = "#457B9D")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  publication_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        axis.title.y = element_text(size = 10, margin = margin(r = 8))) +
  labs(y = "Number of significant DEGs", x = "")

ggsave(file.path(fig_pub_dir, "Figure4_central_finding.pdf"), p_fig4,
       width = 10, height = 6, device = cairo_pdf)
cat("Figure4 fixed (cleaner y-axis label)\n")





