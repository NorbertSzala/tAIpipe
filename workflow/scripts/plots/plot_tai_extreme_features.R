suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

source(snakemake@input[["plotting_utils"]])
source(snakemake@input[["table_validation_utils"]])
source(snakemake@input[["plot_data_utils"]])
source(snakemake@input[["label_utils"]])

plot_cfg <- read_plot_config(snakemake@input[["plotting_config"]])

binary_data <- read_tsv_checked(
  snakemake@input[["binary_summary"]],
  required_columns = c("sample", "tai_group", "feature", "proportion"),
  table_name = "binary_feature_summary.tsv"
)

continuous_data <- read_tsv_checked(
  snakemake@input[["continuous_summary"]],
  required_columns = c("sample", "tai_group", "feature", "median"),
  table_name = "continuous_feature_summary.tsv"
)

group_levels <- c(
  "All genes",
  "Bottom 10%",
  "Top 10%",
  "Bottom 1%",
  "Top 1%"
)

binary_data <- binary_data %>%
  mutate(tai_group = factor(tai_group, levels = group_levels))

continuous_data <- continuous_data %>%
  mutate(tai_group = factor(tai_group, levels = group_levels))

binary_plot <- ggplot(
  binary_data,
  aes(
    x = tai_group,
    y = 100 * proportion,
    fill = tai_group
  )
) +
  geom_boxplot(
    width = plot_cfg$geoms$boxplot$width,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(colour = tai_group),
    width = plot_cfg$geoms$jitter$width,
    size = plot_cfg$geoms$jitter$size,
    alpha = plot_cfg$geoms$jitter$alpha
  ) +
  facet_wrap(vars(feature), ncol = 2) +
  scale_fill_project(
    plot_cfg,
    mapping = "tai_extreme_group",
    levels = group_levels
  ) +
  scale_colour_project(
    plot_cfg,
    mapping = "tai_extreme_group",
    levels = group_levels
  ) +
  labs(
    x = NULL,
    y = "Genes with feature (%)",
    title = "Feature prevalence across tAI groups"
  ) +
  project_theme(plot_cfg) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

continuous_plot <- ggplot(
  continuous_data,
  aes(
    x = tai_group,
    y = median,
    fill = tai_group
  )
) +
  geom_boxplot(
    width = plot_cfg$geoms$boxplot$width,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(colour = tai_group),
    width = plot_cfg$geoms$jitter$width,
    size = plot_cfg$geoms$jitter$size,
    alpha = plot_cfg$geoms$jitter$alpha
  ) +
  facet_wrap(vars(feature), ncol = 2, scales = "free_y") +
  scale_fill_project(
    plot_cfg,
    mapping = "tai_extreme_group",
    levels = group_levels
  ) +
  scale_colour_project(
    plot_cfg,
    mapping = "tai_extreme_group",
    levels = group_levels
  ) +
  labs(
    x = NULL,
    y = "Median value per genome",
    title = "Continuous features across tAI groups"
  ) +
  project_theme(plot_cfg) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

save_plot(binary_plot, snakemake@output[["binary_png"]], plot_cfg, size = "wide")
save_plot(binary_plot, snakemake@output[["binary_pdf"]], plot_cfg, size = "wide")
save_plot(continuous_plot, snakemake@output[["continuous_png"]], plot_cfg, size = "wide")
save_plot(continuous_plot, snakemake@output[["continuous_pdf"]], plot_cfg, size = "wide")