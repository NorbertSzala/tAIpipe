#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# plot_codon_trna_suggestion_set.R
# -----------------------------------------------------------------------------
# Purpose
#   Recreate codon/tRNA plots from script_suggestions:
#   - repair_tRNAscan_other_programmes/tRNAscan-SE_related/barplot_tRNA_anticodons_count.R
#   - cai_enc/scripts/cys_preferences_across_all_fungi.py
#
# Data contract
#   Input: canonical results/tables/codon_profiles.tsv.
#   Required core columns: sample, codon.
#   The script chooses the most appropriate available measure:
#     tRNA copy/weight plots: trna_copy_number > tRNA_weight > trna_absolute_weight
#     codon usage plots: codon_frequency > codon_count > genome_RSCU
#
# Safety
#   If a requested biological signal cannot be reconstructed from the canonical
#   table, no substitute plot is generated. The missing plot is logged as XXXXX.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    codon_profiles = "results/tables/codon_profiles.tsv",
    output_dir = "results/plots/script_suggestions/codon_trna",
    formats = "png,pdf"
  )
  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--")) stop("Unexpected argument: ", args[[i]])
    if (i == length(args)) stop("Missing value for argument: ", args[[i]])
    out[[gsub("-", "_", sub("^--", "", args[[i]]))]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

if (exists("snakemake")) {
  args <- list(
    codon_profiles = snakemake@input[["codon_profiles"]],
    output_dir = snakemake@params[["output_dir"]] %||% "results/plots/script_suggestions/codon_trna",
    formats = paste(snakemake@params[["formats"]] %||% c("png", "pdf"), collapse = ",")
  )
} else {
  args <- parse_args()
}

output_dir <- args$output_dir
formats <- str_split(args$formats, ",", simplify = FALSE)[[1]] |> trimws()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

missing_plots <- tibble(plot = character(), reason = character())
add_missing <- function(plot, reason) {
  missing_plots <<- bind_rows(missing_plots, tibble(plot = plot, reason = reason))
  message("XXXXX skipped ", plot, ": ", reason)
}

save_plot_multi <- function(plot, stem, width = 9, height = 6) {
  for (fmt in formats) {
    ggsave(file.path(output_dir, paste0(stem, ".", fmt)), plot = plot, width = width, height = height, dpi = 300, limitsize = FALSE)
  }
}

pick_col <- function(data, candidates) {
  hit <- intersect(candidates, names(data))
  if (length(hit) == 0L) NA_character_ else hit[[1]]
}

if (!file.exists(args$codon_profiles)) stop("Missing codon_profiles table: ", args$codon_profiles)
codons <- readr::read_tsv(args$codon_profiles, show_col_types = FALSE, progress = FALSE)
if (nrow(codons) == 0L) stop("codon_profiles table is empty")
for (col in c("sample", "codon")) if (!col %in% names(codons)) stop("codon_profiles lacks column: ", col)

# Useful metadata for ordering/faceting. Missing metadata is acceptable.
if (!"amino_acid" %in% names(codons)) codons$amino_acid <- NA_character_
if (!"phylum" %in% names(codons)) codons$phylum <- NA_character_
if (!"species" %in% names(codons)) codons$species <- codons$sample

trna_measure <- pick_col(codons, c("trna_copy_number", "tRNA_weight", "trna_absolute_weight"))
usage_measure <- pick_col(codons, c("codon_frequency", "codon_count", "genome_RSCU"))

# ---- tRNA/codon profile per genome ------------------------------------------

if (!is.na(trna_measure)) {
  trna_df <- codons |>
    mutate(value = suppressWarnings(as.numeric(.data[[trna_measure]]))) |>
    filter(is.finite(value), !is.na(codon))

  if (nrow(trna_df) > 0L) {
    p_heat <- ggplot(trna_df, aes(x = codon, y = sample, fill = value)) +
      geom_tile() +
      facet_grid(rows = vars(phylum), scales = "free_y", space = "free_y") +
      labs(x = "Codon", y = "Genome", fill = trna_measure, title = "tRNA/codon profile per genome") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text.y = element_text(angle = 0))
    save_plot_multi(p_heat, "codon_profiles_per_genome", width = 12, height = 10)

    trna_norm <- trna_df |>
      group_by(sample) |>
      mutate(value_norm = ifelse(max(value, na.rm = TRUE) > 0, value / max(value, na.rm = TRUE), NA_real_)) |>
      ungroup()

    p_heat_norm <- ggplot(trna_norm, aes(x = codon, y = sample, fill = value_norm)) +
      geom_tile() +
      facet_grid(rows = vars(phylum), scales = "free_y", space = "free_y") +
      labs(x = "Codon", y = "Genome", fill = paste0(trna_measure, " normalized"), title = "Normalized tRNA/codon profile per genome") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text.y = element_text(angle = 0))
    save_plot_multi(p_heat_norm, "codon_profiles_per_genome_normalized", width = 12, height = 10)

    variation <- trna_df |>
      group_by(codon, amino_acid) |>
      summarise(
        n = sum(is.finite(value)),
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        cv = ifelse(is.finite(mean) && mean != 0, sd / abs(mean), NA_real_),
        .groups = "drop"
      ) |>
      arrange(desc(cv))

    p_var <- variation |>
      slice_head(n = 40) |>
      ggplot(aes(x = reorder(codon, cv), y = cv, fill = amino_acid)) +
      geom_col() +
      coord_flip() +
      labs(x = NULL, y = "CV across genomes", fill = "Amino acid", title = "Codon/tRNA variation across organisms") +
      theme_minimal(base_size = 10)
    save_plot_multi(p_var, "codon_variation_across_organisms", width = 8, height = 8)

    norm_var <- trna_norm |>
      group_by(codon, amino_acid) |>
      summarise(
        n = sum(is.finite(value_norm)),
        mean = mean(value_norm, na.rm = TRUE),
        sd = sd(value_norm, na.rm = TRUE),
        cv = ifelse(is.finite(mean) && mean != 0, sd / abs(mean), NA_real_),
        .groups = "drop"
      ) |>
      arrange(desc(cv))

    p_var_norm <- norm_var |>
      slice_head(n = 40) |>
      ggplot(aes(x = reorder(codon, cv), y = cv, fill = amino_acid)) +
      geom_col() +
      coord_flip() +
      labs(x = NULL, y = "CV across genomes", fill = "Amino acid", title = "Normalized codon/tRNA variation across organisms") +
      theme_minimal(base_size = 10)
    save_plot_multi(p_var_norm, "codon_variation_across_organisms_normalized", width = 8, height = 8)
  } else {
    add_missing("codon_profiles_per_genome / variation", paste("column has no finite values:", trna_measure))
  }
} else {
  add_missing("codon_profiles_per_genome / variation", "missing tRNA measure: expected trna_copy_number, tRNA_weight, or trna_absolute_weight")
}

# ---- Cys preferences across fungi -------------------------------------------
# Cys DNA codons under the standard nuclear code are TGT and TGC. This plot is
# generated only if codon usage can be represented directly from codon_profiles.

if (!is.na(usage_measure)) {
  cys <- codons |>
    filter(codon %in% c("TGT", "TGC")) |>
    mutate(value = suppressWarnings(as.numeric(.data[[usage_measure]]))) |>
    filter(is.finite(value))

  if (nrow(cys) >= 2L && length(unique(cys$codon)) == 2L) {
    # If counts are used, convert to within-sample Cys codon share. If frequency
    # or RSCU is used, retain original scale but still show per-codon values.
    cys_plot_data <- cys |>
      group_by(sample) |>
      mutate(
        cys_total = sum(value, na.rm = TRUE),
        cys_share = ifelse(cys_total > 0, value / cys_total, NA_real_)
      ) |>
      ungroup()

    y_col <- if (usage_measure == "codon_count") "cys_share" else "value"
    y_label <- if (usage_measure == "codon_count") "Cys codon share" else usage_measure

    p_species <- ggplot(cys_plot_data, aes(x = sample, y = .data[[y_col]], fill = codon)) +
      geom_col(position = "dodge") +
      facet_grid(rows = vars(phylum), scales = "free_y", space = "free_y") +
      labs(x = "Genome", y = y_label, fill = "Cys codon", title = "Cys codon usage per genome") +
      theme_minimal(base_size = 8) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text.y = element_text(angle = 0))
    save_plot_multi(p_species, "codon_usage_per_species", width = 13, height = 10)

    p_distribution <- ggplot(cys_plot_data, aes(x = codon, y = .data[[y_col]], fill = codon)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8) +
      geom_jitter(width = 0.15, size = 0.7, alpha = 0.35) +
      labs(x = "Cys codon", y = y_label, title = "Distribution of Cys codon usage across fungi") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "none")
    save_plot_multi(p_distribution, "codon_usage_distribution", width = 6.5, height = 5.2)
  } else {
    add_missing("Cys codon usage plots", "TGT/TGC rows are missing or have no finite usage values")
  }
} else {
  add_missing("Cys codon usage plots", "missing usage measure: expected codon_frequency, codon_count, or genome_RSCU")
}

if (nrow(missing_plots) == 0L) missing_plots <- tibble(plot = "none", reason = "all possible codon/tRNA plots generated")
readr::write_tsv(missing_plots, file.path(output_dir, "XXXXX_missing_plots.tsv"))
message("Done. Codon/tRNA plots written to: ", output_dir)
