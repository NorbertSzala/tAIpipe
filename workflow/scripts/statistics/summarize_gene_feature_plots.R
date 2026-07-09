suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
})

source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")

statistics_cfg <- snakemake@params[["statistics"]]

binary_features <- statistics_cfg[["binary_features"]]
continuous_features <- statistics_cfg[["gene_covariates"]]

feature_labels <- setNames(binary_features, binary_features)

genes <- read_tsv_checked(
  snakemake@input[["gene_features"]],
  required_columns = unique(c("sample", "gene_id", "tAI", binary_features, continuous_features)),
  table_name = "gene_features.tsv"
)

binary_map <- setNames(binary_features, vapply(binary_features, identity, character(1)))
continuous_map <- setNames(continuous_features, vapply(continuous_features, identity, character(1)))

binary_summary <- summarize_binary_by_group(
  genes,
  group_cols = "sample",
  feature_map = binary_map
)

continuous_summary <- summarize_numeric_by_group(
  genes,
  group_cols = "sample",
  feature_map = continuous_map
)

distribution_summary <- bind_rows(
  binary_summary |> mutate(feature_type = "binary"),
  continuous_summary |> mutate(feature_type = "continuous")
)

gene_tests <- read_tsv_checked(
  snakemake@input[["gene_tests"]],
  table_name = "gene_feature_tests.tsv"
)

write_tsv_safe(distribution_summary, snakemake@output[["distribution_summary"]])
write_tsv_safe(gene_tests, snakemake@output[["effect_summary"]])
