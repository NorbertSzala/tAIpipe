suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
})

source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")


cfg <- snakemake@params[["tai_extremes"]]

group_column <- cfg[["group_column"]]
id_column <- cfg[["id_column"]]
value_column <- cfg[["value_column"]]
tail_fractions <- as.numeric(unlist(cfg[["tail_fractions"]]))

binary_features <- unlist(cfg[["binary_features"]], use.names = TRUE)
continuous_features <- unlist(cfg[["continuous_features"]], use.names = TRUE)

required_columns <- unique(c(
  group_column,
  id_column,
  value_column,
  unname(binary_features),
  unname(continuous_features)
))

genes <- read_tsv_checked(
  snakemake@input[["gene_features"]],
  required_columns = required_columns,
  table_name = "gene_features.tsv"
)

genes <- genes |>
  filter(!is.na(.data[[value_column]]))

format_fraction_label <- function(fraction) {
  format(100 * fraction, scientific = FALSE, trim = TRUE)
}

select_tails <- function(data, fraction) {
  label <- format_fraction_label(fraction)

  data |>
    group_by(.data[[group_column]]) |>
    arrange(.data[[value_column]], .data[[id_column]], .by_group = TRUE) |>
    mutate(
      .n_genes = n(),
      .tail_size = pmax(1L, floor(.n_genes * fraction)),
      .rank = row_number(),
      tai_group = case_when(
        .rank <= .tail_size ~ paste0("Bottom ", label, "%"),
        .rank > .n_genes - .tail_size ~ paste0("Top ", label, "%"),
        TRUE ~ NA_character_
      ),
      tail_fraction = fraction
    ) |>
    filter(.n_genes >= 2, !is.na(tai_group)) |>
    select(-starts_with(".")) |>
    ungroup()
}

tail_membership <- map_dfr(tail_fractions, ~ select_tails(genes, .x))

all_genes <- genes |>
  mutate(
    tai_group = "All genes",
    tail_fraction = NA_real_
  )

analysis_data <- bind_rows(all_genes, tail_membership)

binary_summary <- summarize_binary_by_group(
  analysis_data,
  group_cols = c(group_column, "tai_group", "tail_fraction"),
  feature_map = binary_features
)

continuous_summary <- summarize_numeric_by_group(
  analysis_data,
  group_cols = c(group_column, "tai_group", "tail_fraction"),
  feature_map = continuous_features
)

tail_membership |>
  select(all_of(c(group_column, id_column, value_column)), tai_group, tail_fraction) |>
  write_tsv_safe(snakemake@output[["membership"]])

write_tsv_safe(binary_summary, snakemake@output[["binary_summary"]])
write_tsv_safe(continuous_summary, snakemake@output[["continuous_summary"]])
