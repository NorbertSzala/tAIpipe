suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

source(snakemake@input[["plotting_utils"]])
source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")
source("workflow/scripts/lib/label_utils.R")
source("workflow/scripts/lib/plot_style_helpers.R")

cfg <- read_plot_config(snakemake@input[["plotting_config"]])
summary <- read_tsv_checked(snakemake@input[["distribution_summary"]], table_name = "gene_feature_distribution_summary.tsv")
effects <- read_tsv_checked(snakemake@input[["effect_summary"]], table_name = "gene_feature_effect_summary.tsv")

summary <- summary %>% mutate(feature = clean_plot_label(feature))
effects <- effects %>% mutate(feature = clean_plot_label(feature))

binary <- summary |> filter(feature_type == "binary")
continuous <- summary |> filter(feature_type == "continuous")

binary_p <- binary %>% group_by(feature) %>% summarise(p = safe_kruskal_p(100 * proportion, sample), .groups = "drop")
continuous_p <- continuous %>% group_by(feature) %>% summarise(p = safe_kruskal_p(median, sample), .groups = "drop")

binary_plot <- binary %>%
  ggplot(aes(x = feature, y = 100 * proportion, fill = feature)) +
  geom_violin(trim = FALSE, alpha = 0.65, na.rm = TRUE) +
  geom_boxplot(width = 0.12, outlier.alpha = 0.25, alpha = 0.85, na.rm = TRUE) +
  labs(x = NULL, y = "Genes with feature [%]", title = "Prevalence of binary gene features across genomes") +
  project_theme(cfg) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

continuous_plot <- continuous %>%
  ggplot(aes(x = feature, y = median, fill = feature)) +
  geom_violin(trim = FALSE, alpha = 0.65, na.rm = TRUE) +
  geom_boxplot(width = 0.12, outlier.alpha = 0.25, alpha = 0.85, na.rm = TRUE) +
  labs(
    x = NULL, y = "Per-genome median",
    title = "Distribution of continuous gene features across genomes",
    subtitle = "This is not a GC3s-vs-protein-length biological test; it summarizes each covariate separately. Use regression plots for feature-vs-tAI relationships."
  ) +
  project_theme(cfg) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

# Effect plot: estimates from the statistical model, not raw group means.
# If CI columns are absent, plot only point estimates.
if (all(c("conf_low", "conf_high") %in% names(effects))) {
  effects_plot <- ggplot(effects, aes(x = reorder(as.character(feature), estimate), y = estimate)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.12) +
    geom_point(size = 2)
} else {
  effects_plot <- ggplot(effects, aes(x = reorder(as.character(feature), estimate), y = estimate)) +
    geom_hline(yintercept = 0, linewidth = 0.3) + geom_point(size = 2)
}
effects_plot <- effects_plot + coord_flip() +
  labs(
    x = NULL,
    y = "Model effect estimate",
    title = "Gene-level feature effects",
    subtitle = "Points show estimated change in the modeled tAI response. Zero means no estimated effect; near-zero points are biologically weak even when p-values are small."
  ) + project_theme(cfg)

save_plot(binary_plot, snakemake@output[["binary_png"]], cfg, size = "wide")
save_plot(binary_plot, snakemake@output[["binary_pdf"]], cfg, size = "wide")
save_plot(continuous_plot, snakemake@output[["continuous_png"]], cfg, size = "wide")
save_plot(continuous_plot, snakemake@output[["continuous_pdf"]], cfg, size = "wide")
save_plot(effects_plot, snakemake@output[["effects_png"]], cfg, size = "wide")
save_plot(effects_plot, snakemake@output[["effects_pdf"]], cfg, size = "wide")
