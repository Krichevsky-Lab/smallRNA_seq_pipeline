#!/usr/bin/env Rscript

################################################################################
# 15_analysis.R
#
# Global differential expression analysis using DESeq2 with
# canonical tRF collapsing, followed by a fixed-layout volcano plot
# and PCA of tRNA / mt-tRNA features.
#
# Outputs:
#   analysis/results/deseq2_results.csv
#   analysis/plots/1h.png        (volcano)
#   analysis/plots/1i.png        (PCA)
#   analysis/results/1i.csv      (PCA variance explained)
#
# Designed for publication-quality figures and reproducibility.
################################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(scales)
})

# ==============================================================================
# Paths and directories
# ==============================================================================

counts_file <- "analysis/rna_counts.csv"
meta_file   <- "input/metadata.csv"

results_dir <- "analysis/results"
plots_dir   <- "analysis/plots"

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Experimental design
# ==============================================================================

condition_col <- "stage"

# ==============================================================================
# Load data
# ==============================================================================

counts <- read.csv(counts_file, check.names = FALSE)
meta   <- read.csv(meta_file)

sample_cols <- intersect(colnames(counts), meta$Run)

# ==============================================================================
# Canonical tRF collapsing
# ==============================================================================

trf_classes <- c(
  "5'-tRF", "3'-tRF", "internal-tRF",
  "Mt-5'-tRF", "Mt-3'-tRF", "Mt-internal-tRF"
)

counts <- counts %>%
  mutate(
    genename = if_else(
      biotype_2 %in% trf_classes,
      sub("(-\\d+)+$", "", genename),
      genename
    )
  ) %>%
  group_by(genename) %>%
  summarise(
    across(-all_of(sample_cols), first),
    across(all_of(sample_cols), sum),
    .groups = "drop"
  )

# ==============================================================================
# DESeq2 analysis
# ==============================================================================

counts_mat <- counts %>%
  select(genename, all_of(sample_cols)) %>%
  column_to_rownames("genename") %>%
  as.matrix()

meta <- meta %>%
  filter(Run %in% sample_cols) %>%
  column_to_rownames("Run")

meta[[condition_col]] <- factor(meta[[condition_col]])

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_mat),
  colData   = meta,
  design    = as.formula(paste0("~ ", condition_col))
)

dds <- DESeq(dds)

# ==============================================================================
# Expression filter: all samples ≥ 10 normalized counts
# ==============================================================================

norm_counts <- counts(dds, normalized = TRUE)
keep <- apply(norm_counts, 1, function(x) all(x >= 10))
dds <- dds[keep, ]

# ==============================================================================
# Results and annotation
# ==============================================================================

res_df <- results(dds) %>%
  as.data.frame() %>%
  rownames_to_column("genename") %>%
  left_join(counts, by = "genename")

write.csv(
  res_df,
  file.path(results_dir, "deseq2_results.csv"),
  row.names = FALSE
)

# ==============================================================================
# Volcano plot (tRF / mt-tRF only)
# ==============================================================================

volcano_plot_1h <- function(df, title) {

  df <- df %>%
    filter(
      biotype_2 %in% c(
        "5'-tRF", "3'-tRF", "internal-tRF",
        "Mt-5'-tRF", "Mt-3'-tRF", "Mt-internal-tRF"
      )
    ) %>%
    transmute(
      feature   = genename,
      baseMean  = baseMean,
      log2FC    = log2FoldChange,
      padj      = padj,
      neglogFDR = -log10(padj),
      log2bm    = log2(baseMean + 1),
      sig       = abs(log2FC) > 0.585 & padj < 0.05,
      ismt      = biotype_2 %in% c("Mt-5'-tRF","Mt-3'-tRF","Mt-internal-tRF")
    )

  labels <- tribble(
    ~feature,                     ~label,                        ~xlab, ~ylab,
    "5'-tRF-tRNA-Gly-GCC", "paste(\"5' Gly\"^GCC)",  2.90, 4.6,
    "5'-tRF-tRNA-Glu-CTC", "paste(\"5' Glu\"^CTC)",  2.90, 2.6,
    "5'-tRF-tRNA-Gly-CCC", "paste(\"5' Gly\"^CCC)",  2.90, 7.6,
    "5'-tRF-tRNA-Glu-TTC", "paste(\"5' Glu\"^TTC)",  2.90, 13.6
  )

  df_lab <- inner_join(df, labels, by = "feature")

  ggplot(df, aes(log2FC, neglogFDR)) +

    geom_point(
      aes(size = log2bm),
      shape = 21, fill = "grey75", color = "grey40",
      stroke = 0.3, alpha = 0.9
    ) +

    geom_point(
      data = df %>% filter(sig & log2FC > 0),
      aes(size = log2bm),
      shape = 21, fill = "orange3", color = "orange4",
      stroke = 0.4, alpha = 0.9
    ) +

    geom_point(
      data = df %>% filter(sig & log2FC < 0 & !ismt),
      aes(size = log2bm),
      shape = 21, fill = "#56B4E9", color = "#1F78B4",
      stroke = 0.4, alpha = 0.9
    ) +

    geom_point(
      data = df %>% filter(sig & log2FC < 0 & ismt),
      aes(size = log2bm),
      shape = 21, fill = "#003399", color = "#002266",
      stroke = 0.45, alpha = 0.95
    ) +

    annotate(
      "text", x = -2.5, y = 0.9,
      label = "mt-tRFs",
      color = "#003399", size = 5.6,
      fontface = "bold", hjust = 0
    ) +

    geom_point(
      data = df_lab,
      size = 4, shape = 21,
      fill = "#A50000", color = "black", stroke = 0.6
    ) +

    geom_segment(
      data = df_lab,
      aes(xend = xlab - 0.05, yend = ylab),
      linewidth = 0.5
    ) +

    geom_text(
      data = df_lab,
      aes(x = xlab, y = ylab, label = label),
      parse = TRUE, hjust = 0,
      size = 5.6, fontface = "bold"
    ) +

    geom_vline(xintercept = c(-0.585, 0.585),
               linetype = "dotted", color = "grey40") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dotted", color = "grey40") +

    scale_size_continuous(
      name = "log2(baseMean)",
      breaks = c(8, 12, 16),
      range = c(0.5, 4)
    ) +

    guides(
      size = guide_legend(
        override.aes = list(
          fill   = "grey35",
          color  = "grey35",
          alpha  = 1,
          stroke = 0.4,
          shape  = 21
        )
      )
    ) +

    coord_cartesian(xlim = c(-5, 5), clip = "off") +

    labs(
      title = title,
      x = "log2(Fold Change)",
      y = "-log10(FDR)"
    ) +

    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      panel.grid = element_blank(),
      legend.position = c(0.08, 0.96),
      legend.justification = c(0, 1),
      legend.key = element_blank(),
      plot.margin = margin(10, 120, 10, 25)
    )
}

p_volcano <- volcano_plot_1h(
  res_df,
  "AD vs. Control – Prefrontal Cortex"
)

ggsave(
  file.path(plots_dir, "1h.png"),
  p_volcano,
  width = 7, height = 5, dpi = 300
)

# ==============================================================================
# PCA: tRNA / mt-tRNA
# ==============================================================================

counts_pca <- read.csv(counts_file, check.names = FALSE)
meta_pca   <- read.csv(meta_file)

sample_cols <- intersect(colnames(counts_pca), meta_pca$Run)

counts_trna <- counts_pca %>%
  filter(if_all(all_of(sample_cols), ~ .x >= 10)) %>%
  filter(biotype_1 %in% c("tRNA", "mt-tRNA"))

expr_mat <- counts_trna %>%
  select(all_of(sample_cols)) %>%
  as.matrix()

rownames(expr_mat) <- counts_trna$genename
expr_mat <- log2(expr_mat + 1)

pca <- prcomp(t(expr_mat), scale. = TRUE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("Run") %>%
  mutate(PC1 = -PC1) %>%
  left_join(meta_pca, by = "Run") %>%
  mutate(group = ifelse(stage == "early_stage", "Control", "AD"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, fill = group)) +
  geom_point(shape = 21, size = 3.6, stroke = 1, color = "black") +
  scale_fill_manual(values = c(Control = "#1F77B4", AD = "#D62728")) +
  labs(
    x = paste0("PC1: ", round(var_exp[1] * 100), "% variance"),
    y = paste0("PC2: ", round(var_exp[2] * 100), "% variance")
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave(
  file.path(plots_dir, "1i.png"),
  p_pca,
  width = 6.5, height = 3.5, dpi = 300
)

write.csv(
  tibble(
    PC = paste0("PC", seq_along(var_exp)),
    Variance_Explained = var_exp
  ),
  file.path(results_dir, "pca.csv"),
  row.names = FALSE
)
