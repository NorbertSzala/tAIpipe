#!/usr/bin/env Rscript
# Codon/tRNA profile plots.
#
# Produces:
#   - codon_usage_variability.pdf/png
#   - large tRNA-weight heatmaps with codons on X and organisms on Y.
#
# Codon usage variability uses one explicitly selected usage column. By default
# the default is genome_RSCU. You can still pass --usage-value auto or another
# explicit column to change this behaviour. The selected column is written to
# codon_usage_variability_method.tsv so repeated runs are auditable.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
})

source("workflow/scripts/lib/plot_style_helpers.R")

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    codon_profiles = "results/tables/codon_profiles.tsv",
    output_dir = "results/plots/codon_profiles",
    formats = "png,pdf",
    heatmap_value = "auto",
    usage_value = "genome_RSCU",
    top_n_variable_codons = "30",
    heatmap_large_width = "16.5",
    heatmap_large_height = "23.4",
    heatmap_min_group_n = "5"
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

args <- if (exists("snakemake")) {
  list(
    codon_profiles = snakemake@input[["codon_profiles"]] %||% snakemake@input[[1]],
    output_dir = snakemake@params[["output_dir"]] %||% "results/plots/codon_profiles",
    formats = paste(snakemake@params[["formats"]] %||% c("png", "pdf"), collapse = ","),
    heatmap_value = snakemake@params[["heatmap_value"]] %||% "auto",
    usage_value = snakemake@params[["usage_value"]] %||% "genome_RSCU",
    top_n_variable_codons = as.character(snakemake@params[["top_n_variable_codons"]] %||% 30),
    heatmap_large_width = as.character(snakemake@params[["heatmap_large_width"]] %||% 16.5),
    heatmap_large_height = as.character(snakemake@params[["heatmap_large_height"]] %||% 23.4),
    heatmap_min_group_n = as.character(snakemake@params[["heatmap_min_group_n"]] %||% 5)
  )
} else parse_args()

formats <- strsplit(args$formats, ",", fixed = TRUE)[[1]] |> trimws()
outdir <- args$output_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

top_n_variable_codons <- as.integer(args$top_n_variable_codons)
heatmap_large_width <- as.numeric(args$heatmap_large_width)
heatmap_large_height <- as.numeric(args$heatmap_large_height)
heatmap_min_group_n <- as.integer(args$heatmap_min_group_n)

standard_code <- c(
  TTT="Phe", TTC="Phe", TTA="Leu", TTG="Leu",
  TCT="Ser", TCC="Ser", TCA="Ser", TCG="Ser",
  TAT="Tyr", TAC="Tyr", TAA="Stop", TAG="Stop",
  TGT="Cys", TGC="Cys", TGA="Stop", TGG="Trp",
  CTT="Leu", CTC="Leu", CTA="Leu", CTG="Leu",
  CCT="Pro", CCC="Pro", CCA="Pro", CCG="Pro",
  CAT="His", CAC="His", CAA="Gln", CAG="Gln",
  CGT="Arg", CGC="Arg", CGA="Arg", CGG="Arg",
  ATT="Ile", ATC="Ile", ATA="Ile", ATG="Met",
  ACT="Thr", ACC="Thr", ACA="Thr", ACG="Thr",
  AAT="Asn", AAC="Asn", AAA="Lys", AAG="Lys",
  AGT="Ser", AGC="Ser", AGA="Arg", AGG="Arg",
  GTT="Val", GTC="Val", GTA="Val", GTG="Val",
  GCT="Ala", GCC="Ala", GCA="Ala", GCG="Ala",
  GAT="Asp", GAC="Asp", GAA="Glu", GAG="Glu",
  GGT="Gly", GGC="Gly", GGA="Gly", GGG="Gly"
)

amino_order <- c("Phe", "Leu", "Ile", "Met", "Val", "Ser", "Pro", "Thr", "Ala", "Tyr", "His", "Gln", "Asn", "Lys", "Asp", "Glu", "Cys", "Trp", "Arg", "Gly")

choose_col <- function(data, candidates, requested = "auto", label) {
  if (!identical(requested, "auto")) {
    if (!requested %in% names(data)) stop(label, " column not found: ", requested)
    return(requested)
  }
  hit <- candidates[candidates %in% names(data)][1]
  if (is.na(hit)) stop(label, " column not found. Tried: ", paste(candidates, collapse = ", "))
  hit
}

profiles <- readr::read_tsv(args$codon_profiles, show_col_types = FALSE, progress = FALSE)
required <- c("sample", "codon")
missing <- setdiff(required, names(profiles))
if (length(missing) > 0L) stop("codon_profiles.tsv lacks columns: ", paste(missing, collapse = ", "))

profiles <- profiles %>%
  mutate(
    codon = toupper(as.character(codon)),
    amino_acid = if ("amino_acid" %in% names(.)) as.character(amino_acid) else unname(standard_code[codon]),
    amino_acid = if_else(is.na(amino_acid), "Unknown", amino_acid),
    phylum = if ("phylum" %in% names(.)) clean_label(phylum) else "Unknown phylum",
    lifestyle = if ("lifestyle" %in% names(.)) clean_label(lifestyle) else "Unknown lifestyle"
  ) %>%
  filter(amino_acid != "Stop")

heatmap_col <- choose_col(
  profiles,
  c("tRNA_weight", "trna_weight", "trna_absolute_weight", "trna_copy_number", "relative_trna"),
  args$heatmap_value,
  "Heatmap value"
)
usage_col <- choose_col(
  profiles,
  c("genome_RSCU", "RSCU", "rscu", "codon_frequency", "frequency", "codon_count"),
  args$usage_value,
  "Codon usage variability value"
)

profiles <- profiles %>%
  mutate(
    heatmap_value = suppressWarnings(as.numeric(.data[[heatmap_col]])),
    usage_value = suppressWarnings(as.numeric(.data[[usage_col]]))
  )

readr::write_tsv(
  tibble(
    heatmap_value_column = heatmap_col,
    codon_usage_variability_column = usage_col,
    interpretation = ifelse(grepl("RSCU|rscu", usage_col),
      "CV of RSCU across genomes: variability of synonymous codon preference.",
      "CV of codon frequency/count across genomes: variability of raw codon abundance, more composition-dependent."
    )
  ),
  file.path(outdir, "codon_usage_variability_method.tsv")
)

# Main codon-usage variability figure.
cv_tbl <- profiles %>%
  filter(is.finite(usage_value)) %>%
  group_by(codon, amino_acid) %>%
  summarise(
    mean_usage = mean(usage_value, na.rm = TRUE),
    sd_usage = sd(usage_value, na.rm = TRUE),
    cv = sd_usage / abs(mean_usage),
    n_genomes = n_distinct(sample),
    .groups = "drop"
  ) %>%
  filter(is.finite(cv), n_genomes >= 3L) %>%
  arrange(desc(cv)) %>%
  mutate(codon_label = paste0(codon, " (", amino_acid, ")"))

readr::write_tsv(cv_tbl, file.path(outdir, "codon_usage_variability.tsv"))

p_cv <- cv_tbl %>%
  slice_head(n = min(top_n_variable_codons, nrow(cv_tbl))) %>%
  mutate(codon_label = forcats::fct_reorder(codon_label, cv)) %>%
  ggplot(aes(x = codon_label, y = cv)) +
  geom_col(width = 0.75, fill = "#666666") +
  coord_flip() +
  labs(
    x = NULL,
    y = paste0("Coefficient of variation of ", usage_col),
    title = "Most variable codons across genomes",
    subtitle = wrap_text(paste0("CV = SD / |mean| across genomes. This plot uses exactly one selected column: ", usage_col, "."), 90)
  ) +
  theme_minimal(base_size = 11)

save_plot_pair(p_cv, "codon_usage_variability", outdir, 7.8, 6.4, formats)

prepare_heatmap_df <- function(group_col, min_count = NULL) {
  group_col <- rlang::as_string(rlang::ensym(group_col))
  df <- profiles %>%
    filter(is.finite(heatmap_value)) %>%
    mutate(
      group_raw = .data[[group_col]],
      sample_short = unname(short_sample_labels(sample)),
      amino_acid = factor(amino_acid, levels = amino_order),
      codon_label = paste0(codon, "\n", amino_acid)
    )

  if (!is.null(min_count)) {
    tab <- df %>% distinct(sample, group_raw) %>% count(group_raw, name = "n")
    keep <- tab %>% filter(n >= min_count) %>% pull(group_raw)
    other_members <- tab %>% filter(!group_raw %in% keep) %>% pull(group_raw)
    df <- df %>% mutate(group_raw = if_else(group_raw %in% keep, group_raw, "Other"))
    attr(df, "other_members") <- other_members
  }

  group_labels <- df %>% distinct(sample, group_raw) %>% count(group_raw, name = "n") %>%
    mutate(group_label = paste0(clean_label(group_raw), " (", n, ")"))
  label_map <- stats::setNames(group_labels$group_label, group_labels$group_raw)

  codon_levels <- df %>%
    distinct(amino_acid, codon, codon_label) %>%
    arrange(amino_acid, codon) %>%
    pull(codon_label)

  sample_levels <- df %>%
    distinct(group_raw, sample_short) %>%
    arrange(group_raw, sample_short) %>%
    pull(sample_short)

  df %>%
    mutate(
      group_label = label_map[group_raw],
      group_label = stringr::str_wrap(group_label, width = 24),
      sample_short = factor(sample_short, levels = unique(sample_levels)),
      codon_label = factor(codon_label, levels = unique(codon_levels)),
      amino_acid = factor(amino_acid, levels = amino_order)
    )
}

make_large_heatmap <- function(group_col, suffix, min_count = NULL) {
  group_name <- rlang::as_string(rlang::ensym(group_col))
  df <- prepare_heatmap_df({{ group_col }}, min_count = min_count)

  other <- attr(df, "other_members")
  subtitle <- if (!is.null(other) && length(other) > 0L) {
    wrap_text(paste0("Other contains: ", paste(clean_label(other), collapse = ", "), "."), 110)
  } else {
    wrap_text("Samples are ordered within each group. Codons are grouped by encoded amino acid. Full sample order is saved as TSV.", 110)
  }

  readr::write_tsv(
    df %>% distinct(sample, sample_short, group_label) %>% arrange(group_label, sample_short),
    file.path(outdir, paste0("sample_order_by_", group_name, suffix, ".tsv"))
  )

  p <- ggplot(df, aes(x = codon_label, y = sample_short, fill = heatmap_value)) +
    geom_tile(width = 0.97, height = 0.97) +
    facet_grid(group_label ~ amino_acid, scales = "free", space = "free", switch = "y") +
    scale_fill_viridis_c(option = "C", begin = 0.12, end = 0.88, na.value = "#D9D9D9") +
    labs(
      x = "Codon grouped by encoded amino acid",
      y = "Genome",
      fill = heatmap_col,
      title = paste0("tRNA weights by codon and genome, grouped by ", clean_label(group_name)),
      subtitle = subtitle
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5),
      axis.text.y = element_text(size = 8),
      strip.text.x = element_text(size = 8.5, face = "bold"),
      strip.text.y.left = element_text(size = 8.5, face = "bold", angle = 0, hjust = 1),
      strip.placement = "outside",
      legend.position = "right",
      plot.subtitle = element_text(lineheight = 0.95),
      panel.spacing.x = grid::unit(0.08, "lines"),
      panel.spacing.y = grid::unit(0.20, "lines")
    )

  stem <- paste0("trna_weights_heatmap_by_", group_name, suffix, "_large_codon_x")
  save_plot_pair(p, stem, outdir, heatmap_large_width, heatmap_large_height, formats)
}

make_large_heatmap(phylum, "")
make_large_heatmap(lifestyle, "")
make_large_heatmap(phylum, paste0("_count", heatmap_min_group_n), min_count = heatmap_min_group_n)
make_large_heatmap(lifestyle, paste0("_count", heatmap_min_group_n), min_count = heatmap_min_group_n)

message("Done: ", outdir)
