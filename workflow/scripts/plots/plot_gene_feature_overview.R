# Creates the canonical gene-feature figures from per-genome distribution and
# effect summaries, including split violins, correlation summaries and effects.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

source(snakemake@input[["plotting_utils"]])
source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")
source("workflow/scripts/lib/label_utils.R")
source("workflow/scripts/lib/plot_style_helpers.R")

cfg <- read_plot_config(snakemake@input[["plotting_config"]])
summary <- read_tsv_checked(
  snakemake@input[["distribution_summary"]],
  required_columns = c("feature_type", "sample", "feature", "value"),
  table_name = "gene_feature_distribution_summary.tsv"
)
effects <- read_tsv_checked(snakemake@input[["effect_summary"]], table_name = "gene_feature_effect_summary.tsv")

summary <- summary %>%
  mutate(
    feature = clean_plot_label(feature),
    value = suppressWarnings(as.numeric(value))
  )
effects <- effects %>% mutate(feature = clean_plot_label(feature))

binary <- summary %>%
  filter(feature_type == "binary", is.finite(value), !is.na(status)) %>%
  mutate(status = factor(status, levels = c("Absent", "Present")))
continuous <- summary %>% filter(feature_type == "continuous", is.finite(value))

# Split violins compare presence/absence within the same biological feature;
# unrelated annotations are separated into facets rather than overlaid.
if (nrow(binary) > 0L) {
  counts_long <- binary %>%
    group_by(feature, status) %>%
    summarise(n_genes = sum(suppressWarnings(as.numeric(n_genes)), na.rm = TRUE), .groups = "drop") %>%
    mutate(
      facet_label = feature,
      label_x = if_else(as.character(status) == "Absent", 0.72, 1.28),
      count_label = paste0(
        as.character(status), "\nn = ",
        format(n_genes, scientific = FALSE, trim = TRUE, big.mark = ",")
      )
    )
  binary <- binary %>%
    mutate(
      facet_label = feature,
      split_side = split_side_from_level(status, "Absent")
    )

  binary_plot <- ggplot(binary, aes(x = 1, y = value, fill = status,
                                    group = interaction(feature, status), split_side = split_side)) +
    geom_split_violin_project(trim = TRUE, alpha = 0.74, width = 0.94,
                              scale = "area", colour = "grey30", linewidth = 0.22) +
    add_split_violin_boxplots() +
    geom_text(
      data = counts_long,
      aes(x = label_x, y = -Inf, label = count_label),
      inherit.aes = FALSE,
      vjust = 2.30,
      size = 4.05,
      lineheight = 1.05,
      colour = "grey25"
    ) +
    facet_grid(. ~ facet_label, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = binary_grey_values(c("Absent", "Present")), drop = FALSE) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    coord_cartesian(ylim = c(0, NA), clip = "off") +
    labs(
      x = NULL,
      y = "Per-genome median tAI",
      fill = NULL,
      title = "tAI by binary gene features",
      subtitle = "Each violin contains one median per genome and status; separate labels below each half report total gene counts."
    ) +
    project_theme(cfg) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 13.5),
      strip.text = element_text(face = "bold", size = 13.2),
      legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 12.8),
      plot.subtitle = element_text(size = 12.5),
      plot.margin = margin(8, 9, 46, 9)
    )
} else {
  binary_plot <- ggplot() + theme_void() + labs(title = "No binary-feature tAI summaries available")
}

# Per-genome correlations have a common, interpretable [-1, 1] scale. The two
# configured covariates use neutral grey/navy colours for high contrast.
if (nrow(continuous) > 0L) {
  features <- unique(as.character(continuous$feature))
  neutral_values <- setNames(rep(c("#6F6F6F", "#1F4E79"), length.out = length(features)), features)
  continuous_plot <- continuous %>%
    mutate(feature = factor(feature, levels = features)) %>%
    ggplot(aes(x = feature, y = value, fill = feature)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "grey45") +
    geom_violin(trim = TRUE, alpha = 0.74, width = 0.92, scale = "width") +
    geom_boxplot(width = 0.10, outlier.alpha = 0.25, alpha = 0.90, fill = "white") +
    scale_fill_manual(values = neutral_values, guide = "none") +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    labs(
      x = NULL,
      y = "Within-genome Spearman rho (tAI vs feature)",
      title = "tAI associations with continuous gene features",
      subtitle = "Each observation is one genome; the dashed line denotes no monotonic association."
    ) +
    project_theme(cfg) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          plot.title = element_text(size = 15, face = "bold"))
} else {
  continuous_plot <- ggplot() + theme_void() + labs(title = "No continuous-feature correlation summaries available")
}

# Effects: per genome, median tAI_z(present) minus median tAI_z(absent). Point and
# interval summarize the distribution across genomes with a signed-rank model.
effects <- effects %>%
  mutate(
    estimate = suppressWarnings(as.numeric(estimate)),
    conf_low = suppressWarnings(as.numeric(conf_low)),
    conf_high = suppressWarnings(as.numeric(conf_high)),
    q_value = if ("q_value" %in% names(.)) suppressWarnings(as.numeric(q_value)) else NA_real_
  ) %>%
  filter(is.finite(estimate)) %>%
  arrange(estimate) %>%
  mutate(feature = factor(feature, levels = unique(feature)))

if (all(c("conf_low", "conf_high") %in% names(effects))) {
  effects_plot <- ggplot(effects, aes(x = feature, y = estimate)) +
    geom_hline(yintercept = 0, linewidth = 0.45, linetype = "dashed", colour = "grey45") +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.12,
                  linewidth = 0.75, colour = "#1F4E79") +
    geom_point(size = 3.0, shape = 21, fill = "white", colour = "#1F4E79", stroke = 1.0)
} else {
  effects_plot <- ggplot(effects, aes(x = feature, y = estimate)) +
    geom_hline(yintercept = 0, linewidth = 0.45, linetype = "dashed", colour = "grey45") +
    geom_point(size = 3.0, colour = "#1F4E79")
}

effects_plot <- effects_plot +
  coord_flip() +
  labs(
    x = NULL,
    y = "Difference in median tAI z-score [present - absent]",
    title = "Genome-level effects of binary gene features",
    subtitle = "Point: Hodges-Lehmann pseudomedian; interval: 95% Wilcoxon CI; dashed line: no common shift"
  ) +
  project_theme(cfg) +
  theme(plot.subtitle = element_text(size = 10), plot.title = element_text(size = 15, face = "bold"))

save_plot(binary_plot, snakemake@output[["binary_png"]], cfg, size = "go_wide_tall")
save_plot(binary_plot, snakemake@output[["binary_pdf"]], cfg, size = "go_wide_tall")
save_plot(continuous_plot, snakemake@output[["continuous_png"]], cfg, size = "wide")
save_plot(continuous_plot, snakemake@output[["continuous_pdf"]], cfg, size = "wide")
save_plot(effects_plot, snakemake@output[["effects_png"]], cfg, size = "wide")
save_plot(effects_plot, snakemake@output[["effects_pdf"]], cfg, size = "wide")
