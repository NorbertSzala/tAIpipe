suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")


codons <- read_tsv_checked(
  snakemake@input[["codon_profiles"]],
  required_columns = c("sample"),
  table_name = "codon_profiles.tsv"
)

# This script is intentionally conservative because codon profile schemas vary.
# It preserves available long-format codon data if columns exist.
possible_measure_cols <- intersect(
  c("tRNA_weight", "relative_trna", "RSCU", "rscu", "codon_count", "frequency"),
  names(codons)
)

if (length(possible_measure_cols) == 0) {
  stop(
    "codon_profiles.tsv does not contain a recognized measure column. ",
    "Expected one of: tRNA_weight, relative_trna, RSCU, rscu, codon_count, frequency.",
    call. = FALSE
  )
}

measure_col <- possible_measure_cols[[1]]

heatmap_cols <- intersect(c("sample", "amino_acid", "codon", measure_col), names(codons))
require_columns(codons, heatmap_cols, "codon_profiles.tsv")

heatmap <- codons |>
  select(all_of(heatmap_cols)) |>
  rename(value = all_of(measure_col))

variability <- heatmap |>
  group_by(dplyr::across(any_of(c("amino_acid", "codon")))) |>
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    cv = ifelse(isTRUE(all.equal(mean, 0)), NA_real_, sd / abs(mean)),
    .groups = "drop"
  )

# Placeholder for full-vs-reference comparison if the input table already contains profile_type.
if ("profile_type" %in% names(codons)) {
  reference_comparison <- codons |>
    select(any_of(c("sample", "amino_acid", "codon", "profile_type", measure_col))) |>
    rename(value = all_of(measure_col))
} else {
  reference_comparison <- tibble::tibble()
}

write_tsv_safe(heatmap, snakemake@output[["heatmap"]])
write_tsv_safe(variability, snakemake@output[["variability"]])
write_tsv_safe(reference_comparison, snakemake@output[["reference_comparison"]])
