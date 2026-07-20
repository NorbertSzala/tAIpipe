suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")


top_n <- as.integer(snakemake@params[["top_n"]])

go <- read_tsv_checked(
  snakemake@input[["go_enrichment"]],
  table_name = "go_enrichment.tsv"
)

p_col <- intersect(c("padj", "fdr", "q_value", "adjusted_p_value"), names(go))[[1]]
if (is.na(p_col)) {
  stop("No adjusted p-value column found in GO enrichment table.", call. = FALSE)
}

top_terms <- select_top_terms(go, p_col = p_col, n = top_n)

if ("namespace" %in% names(go)) {
  namespace_summary <- go |>
    group_by(namespace) |>
    summarise(
      n_terms = n(),
      n_significant = sum(.data[[p_col]] < 0.05, na.rm = TRUE),
      min_padj = min(.data[[p_col]], na.rm = TRUE),
      .groups = "drop"
    )
} else {
  namespace_summary <- tibble::tibble()
}

write_tsv_safe(top_terms, snakemake@output[["top_terms"]])
write_tsv_safe(namespace_summary, snakemake@output[["namespace_summary"]])
