suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")


statistics_cfg <- snakemake@params[["statistics"]]
metrics <- statistics_cfg[["genome_metrics"]]
groups <- statistics_cfg[["group_variables"]]

genomes <- read_tsv_checked(
  snakemake@input[["genome_summary"]],
  required_columns = unique(c("sample", metrics, groups)),
  table_name = "genome_summary.tsv"
)

metric_summary <- genomes |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(metrics),
    names_to = "metric",
    values_to = "value"
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(groups),
    names_to = "group_variable",
    values_to = "group_value"
  ) |>
  group_by(group_variable, group_value, metric) |>
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    q1 = as.numeric(quantile(value, 0.25, na.rm = TRUE)),
    q3 = as.numeric(quantile(value, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )

genome_tests <- read_tsv_checked(
  snakemake@input[["genome_tests"]],
  table_name = "genome_group_tests.tsv"
)

write_tsv_safe(metric_summary, snakemake@output[["metric_summary"]])
write_tsv_safe(genome_tests, snakemake@output[["effect_summary"]])
