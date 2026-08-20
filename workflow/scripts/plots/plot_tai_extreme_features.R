suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

source(snakemake@input[["plotting_utils"]])
source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")
source("workflow/scripts/lib/plot_style_helpers.R")

cfg <- read_plot_config(snakemake@input[["plotting_config"]])

binary <- read_tsv_checked(snakemake@input[["binary_summary"]], table_name = "binary_feature_summary.tsv")
continuous <- read_tsv_checked(snakemake@input[["continuous_summary"]], table_name = "continuous_feature_summary.tsv")

format_tail_group <- function(x) {
  x <- as.character(x)
  x <- gsub("Bottom ", "Bottom ", x)
  x <- gsub("Top ", "Top ", x)
  factor(x, levels = c("All genes", "Bottom 10%", "Top 10%", "Bottom 1%", "Top 1%"))
}

binary <- binary %>%
  mutate(
    feature = clean_plot_label(feature),
    tai_group = format_tail_group(tai_group),
    proportion_pct = 100 * as.numeric(proportion)
  )
continuous <- continuous %>%
  mutate(
    feature = clean_plot_label(feature),
    tai_group = format_tail_group(tai_group),
    median = as.numeric(median)
  )

# One p-value per feature: Kruskal-Wallis across per-genome tail-group summaries.
binary_p <- binary %>% group_by(feature) %>% summarise(p = safe_kruskal_p(proportion_pct, tai_group), .groups = "drop")
continuous_p <- continuous %>% group_by(feature) %>% summarise(p = safe_kruskal_p(median, tai_group), .groups = "drop")

binary_plot <- binary %>% left_join(binary_p, by = "feature") %>%
  ggplot(aes(x = tai_group, y = proportion_pct, fill = tai_group)) +
  geom_violin(trim = FALSE, alpha = 0.65, na.rm = TRUE) +
  geom_boxplot(width = 0.10, outlier.alpha = 0.25, alpha = 0.85, na.rm = TRUE) +
  geom_text(
    data = binary_p, aes(x = Inf, y = Inf, label = sig_label(p)),
    inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 3.6
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  facet_wrap(vars(feature), scales = "fixed") +
  labs(
    x = NULL,
    y = "Genes with feature [%]",
    title = tex_label("Binary features in $\\mathrm{tAI}$ distribution tails")
  ) +
  project_theme(cfg) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 13.5),
    axis.text.y = element_text(size = 13.5),
    axis.title = element_text(size = 14.5),
    strip.text = element_text(size = 14.5, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "none"
  )

continuous_plot <- continuous %>% left_join(continuous_p, by = "feature") %>%
  ggplot(aes(x = tai_group, y = median, fill = tai_group)) +
  geom_violin(trim = FALSE, alpha = 0.65, na.rm = TRUE) +
  geom_boxplot(width = 0.10, outlier.alpha = 0.25, alpha = 0.85, na.rm = TRUE) +
  geom_text(
    data = continuous_p, aes(x = Inf, y = Inf, label = sig_label(p)),
    inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 3.6
  ) +
  facet_wrap(vars(feature), scales = "free_y") +
  labs(
    x = NULL,
    y = "Per-genome median feature value",
    title = tex_label("Continuous features in $\\mathrm{tAI}$ distribution tails")
  ) +
  project_theme(cfg) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 13.5),
    axis.text.y = element_text(size = 13.5),
    axis.title = element_text(size = 14.5),
    strip.text = element_text(size = 14.5, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "none"
  )

save_plot(binary_plot, snakemake@output[["binary_png"]], cfg, size = "wide")
save_plot(binary_plot, snakemake@output[["binary_pdf"]], cfg, size = "wide")
save_plot(continuous_plot, snakemake@output[["continuous_png"]], cfg, size = "wide")
save_plot(continuous_plot, snakemake@output[["continuous_pdf"]], cfg, size = "wide")
