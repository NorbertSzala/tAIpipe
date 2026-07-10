suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")
source("workflow/scripts/lib/plot_style_helpers.R")

codons <- read_tsv_checked(
  snakemake@input[["codon_profiles"]],
  required_columns = c("sample", "codon"),
  table_name = "codon_profiles.tsv"
)

# Normalize optional metadata columns so downstream plotting code does not fail
# when an older codon_profiles.tsv lacks one of them.
if (!"species" %in% names(codons)) codons$species <- codons$sample
if (!"phylum" %in% names(codons)) codons$phylum <- "Unknown"
if (!"lifestyle" %in% names(codons)) codons$lifestyle <- "Unknown"
if (!"amino_acid" %in% names(codons)) codons$amino_acid <- codon_to_aa1(codons$codon)

# Separate biological quantities:
#   - heatmap_value: tRNA/codon adaptation weights, preferably tRNA_weight.
#   - usage_value: observed codon usage, preferably codon_frequency.
# This avoids the earlier ambiguity where the first available numeric column could
# be used for both plots.
heatmap_measure <- intersect(
  c("tRNA_weight", "trna_absolute_weight", "relative_trna", "trna_copy_number"),
  names(codons)
)[1]
usage_measure <- intersect(
  c("codon_frequency", "genome_RSCU", "RSCU", "rscu", "codon_count"),
  names(codons)
)[1]

if (is.na(heatmap_measure)) {
  stop("No recognized tRNA/adaptation measure for heatmap. Expected tRNA_weight, trna_absolute_weight, relative_trna, or trna_copy_number.", call. = FALSE)
}
if (is.na(usage_measure)) {
  stop("No recognized codon-usage measure for variability. Expected codon_frequency, genome_RSCU, RSCU, rscu, or codon_count.", call. = FALSE)
}

meta_cols <- intersect(
  c("sample", "species", "accession", "phylum", "lifestyle", "amino_acid", "aa_code", "subfam", "codon", "is_stop_codon"),
  names(codons)
)

heatmap <- codons %>%
  select(all_of(meta_cols), value = all_of(heatmap_measure)) %>%
  mutate(
    value = suppressWarnings(as.numeric(value)),
    amino_acid = dplyr::coalesce(as.character(amino_acid), codon_to_aa1(codon)),
    species_label = short_species_label(species, sample),
    measure = heatmap_measure
  )

usage <- codons %>%
  select(all_of(meta_cols), usage_value = all_of(usage_measure)) %>%
  mutate(
    usage_value = suppressWarnings(as.numeric(usage_value)),
    amino_acid = dplyr::coalesce(as.character(amino_acid), codon_to_aa1(codon)),
    measure = usage_measure
  )

# CV = sd / abs(mean). Values > 1 are possible and indicate that variability
# across genomes exceeds the mean usage. This is common for rare codons with many
# near-zero values; therefore we report the mean as a diagnostic and filter only
# for at least 10 non-missing genomes.
variability <- usage %>%
  group_by(across(any_of(c("amino_acid", "codon")))) %>%
  summarise(
    n = sum(is.finite(usage_value)),
    mean = mean(usage_value, na.rm = TRUE),
    sd = sd(usage_value, na.rm = TRUE),
    cv = ifelse(is.finite(mean) && abs(mean) > 0, sd / abs(mean), NA_real_),
    measure = first(measure),
    .groups = "drop"
  ) %>%
  filter(n >= 10L)

reference_comparison <- if ("profile_type" %in% names(codons)) {
  codons %>%
    select(any_of(c("sample", "species", "phylum", "lifestyle", "amino_acid", "codon", "profile_type", heatmap_measure))) %>%
    rename(value = all_of(heatmap_measure))
} else {
  tibble::tibble()
}

write_tsv_safe(heatmap, snakemake@output[["heatmap"]])
write_tsv_safe(variability, snakemake@output[["variability"]])
write_tsv_safe(reference_comparison, snakemake@output[["reference_comparison"]])
