suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
})

source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")

statistics_cfg <- snakemake@params[["statistics"]]
binary_features <- statistics_cfg[["binary_features"]]
continuous_features <- statistics_cfg[["gene_covariates"]]

required <- unique(c("sample", "gene_id", "tAI", binary_features, continuous_features))
genes <- read_tsv_checked(
  snakemake@input[["gene_features"]],
  required_columns = required,
  table_name = "gene_features.tsv"
) %>%
  mutate(tAI = suppressWarnings(as.numeric(tAI)))

# Binary plot summary: one raw-tAI median for each genome × feature × status.
# This preserves the genome as the independent visual unit and keeps raw tAI on
# its natural non-negative scale.
binary_summary <- purrr::map_dfr(binary_features, function(feature_name) {
  genes %>%
    transmute(
      sample = as.character(sample),
      feature = feature_name,
      status = coerce_binary_feature(.data[[feature_name]], feature_name),
      tAI = tAI
    ) %>%
    filter(!is.na(status), is.finite(tAI)) %>%
    group_by(sample, feature, status) %>%
    summarise(
      n_genes = n(),
      value = median(tAI, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      status = if_else(status, "Present", "Absent"),
      feature_type = "binary",
      statistic = "per_genome_median_tAI"
    )
})

# Continuous plot summary: one within-genome Spearman coefficient for each
# continuous feature. Spearman is robust to nonlinear monotonic relationships
# and has the same value for raw tAI and any monotonic rescaling of tAI.
continuous_summary <- purrr::map_dfr(continuous_features, function(feature_name) {
  genes %>%
    transmute(
      sample = as.character(sample),
      feature = feature_name,
      tAI = tAI,
      x = suppressWarnings(as.numeric(.data[[feature_name]]))
    ) %>%
    filter(is.finite(tAI), is.finite(x)) %>%
    group_by(sample, feature) %>%
    summarise(
      n_genes = n(),
      value = if (n() >= 10L && n_distinct(x) >= 3L && n_distinct(tAI) >= 3L) {
        suppressWarnings(stats::cor(tAI, x, method = "spearman", use = "complete.obs"))
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      status = NA_character_,
      feature_type = "continuous",
      statistic = "within_genome_spearman_rho_tAI_vs_feature"
    )
})

distribution_summary <- bind_rows(binary_summary, continuous_summary) %>%
  select(feature_type, statistic, sample, feature, status, n_genes, value)

gene_tests <- read_tsv_checked(
  snakemake@input[["gene_tests"]],
  table_name = "gene_feature_tests.tsv"
)

write_tsv_safe(distribution_summary, snakemake@output[["distribution_summary"]])
write_tsv_safe(gene_tests, snakemake@output[["effect_summary"]])
