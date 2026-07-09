suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

source(snakemake@input[["plotting_utils"]])
source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")
source("workflow/scripts/lib/label_utils.R")

cfg <- read_plot_config(snakemake@input[["plotting_config"]])

summary <- read_tsv_checked(
  snakemake@input[["distribution_summary"]],
  table_name = "gene_feature_distribution_summary.tsv"
)

effects <- read_tsv_checked(
  snakemake@input[["effect_summary"]],
  table_name = "gene_feature_effect_summary.tsv"
)

binary <- summary |> filter(feature_type == "binary")
continuous <- summary |> filter(feature_type == "continuous")

binary_plot <- ggplot(binary, aes(x = feature, y = 100 * proportion)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.45) +
  labs(x = NULL, y = "Genes with feature (%)", title = "Prevalence of binary gene features") +
  project_theme(cfg) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

continuous_plot <- ggplot(continuous, aes(x = feature, y = median)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.45) +
  labs(x = NULL, y = "Per-genome median", title = "Distribution of continuous gene features") +
  project_theme(cfg) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

effects_plot <- ggplot(effects, aes(x = reorder(as.character(feature), estimate), y = estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point() +
  coord_flip() +
  labs(x = NULL, y = "Effect estimate", title = "Gene-level feature effects") +
  project_theme(cfg)

save_plot(binary_plot, snakemake@output[["binary_png"]], cfg, size = "wide")
save_plot(binary_plot, snakemake@output[["binary_pdf"]], cfg, size = "wide")
save_plot(continuous_plot, snakemake@output[["continuous_png"]], cfg, size = "wide")
save_plot(continuous_plot, snakemake@output[["continuous_pdf"]], cfg, size = "wide")
save_plot(effects_plot, snakemake@output[["effects_png"]], cfg, size = "wide")
save_plot(effects_plot, snakemake@output[["effects_pdf"]], cfg, size = "wide")
