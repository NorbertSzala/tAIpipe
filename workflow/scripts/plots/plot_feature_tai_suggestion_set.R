#!/usr/bin/env Rscript
# Fast, report-oriented feature-vs-tAI plots.
#
# Design choices:
#   - Histograms use all genes but plot gene counts in millions.
#   - Violin plots use bounded stratified sampling for speed only; tests are not
#     computed on the sampled rows.
#   - Tests use per-genome medians as the unit of analysis to avoid inflated
#     significance caused by millions of non-independent gene rows.
#   - Legend/axis labels include total gene counts so proportions are visible.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
})

source("workflow/scripts/lib/plot_style_helpers.R")

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    gene_table = "results/tables/gene_features.tsv",
    output_dir = "results/plots/script_suggestions/features_tai",
    metric = "tAI",
    max_violin_rows_per_group = "80000",
    max_scatter_rows = "100000",
    seed = "1",
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

args <- if (exists("snakemake")) {
  list(
    gene_table = snakemake@input[["gene_features"]],
    output_dir = snakemake@params[["output_dir"]] %||% "results/plots/script_suggestions/features_tai",
    metric = snakemake@params[["metric"]] %||% "tAI",
    max_violin_rows_per_group = as.character(snakemake@params[["max_violin_rows_per_group"]] %||% 80000),
    max_scatter_rows = as.character(snakemake@params[["max_scatter_rows"]] %||% 100000),
    seed = as.character(snakemake@params[["seed"]] %||% 1),
    formats = paste(snakemake@params[["formats"]] %||% c("png", "pdf"), collapse = ",")
  )
} else parse_args()

output_dir <- args$output_dir
formats <- strsplit(args$formats, ",", fixed = TRUE)[[1]] |> trimws()
metric <- args$metric
max_violin_rows_per_group <- as.integer(args$max_violin_rows_per_group)
max_scatter_rows <- as.integer(args$max_scatter_rows)
set.seed(as.integer(args$seed))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

missing_plots <- tibble(plot = character(), reason = character())
add_missing <- function(plot, reason) {
  missing_plots <<- bind_rows(missing_plots, tibble(plot = plot, reason = reason))
  message("XXXXX skipped ", plot, ": ", reason)
}

if (!file.exists(args$gene_table)) stop("Missing gene table: ", args$gene_table)
genes <- readr::read_tsv(args$gene_table, show_col_types = FALSE, progress = FALSE)
if (!metric %in% names(genes)) stop("Metric column not found: ", metric)
if (!"sample" %in% names(genes)) stop("gene_features.tsv lacks sample column")

# Backwards-compatible alias: older plot requests used pfam_lcr_present.
if ("pfam_lcr_overlap_present" %in% names(genes) && !"pfam_lcr_present" %in% names(genes)) {
  genes$pfam_lcr_present <- genes$pfam_lcr_overlap_present
}

genes <- genes %>%
  mutate(.metric = suppressWarnings(as.numeric(.data[[metric]]))) %>%
  filter(is.finite(.metric))

binary_features <- tibble::tribble(
  ~stem, ~label, ~column,
  "domains", "PFAM domains", "pfam_present",
  "LCR_presence", "LCR", "lcr_present",
  "signal_presence", "Signal peptide", "signal_peptide_present",
  "TM_presence", "Transmembrane region", "tm_present",
  "PFAM_LCR", "PFAM-LCR overlap", "pfam_lcr_present"
)

continuous_features <- tibble::tribble(
  ~stem, ~label, ~column,
  "protein_length", "Protein length", "protein_length_aa",
  "LCR_length", "Total LCR length", "lcr_total_length",
  "TM_length", "Total TM length", "tm_total_length",
  "LCR_number", "LCR count", "lcr_count",
  "TM_count", "TM count", "tm_count"
)

sample_per_group <- function(df, group_col, max_n) {
  # dplyr::slice_sample(n = min(n(), max_n)) is invalid because n must be a
  # constant. group_modify computes a constant n separately for every group.
  df %>%
    group_by(.data[[group_col]]) %>%
    group_modify(function(.x, .y) {
      k <- min(nrow(.x), max_n)
      if (k <= 0L) return(.x[0, , drop = FALSE])
      dplyr::slice_sample(.x, n = k, replace = FALSE)
    }) %>%
    ungroup()
}

make_binary_labels <- function(df) {
  tab <- df %>% count(present, name = "n")
  stats::setNames(paste0(as.character(tab$present), " (n=", format(tab$n, big.mark = ","), ")"), tab$present)
}

plot_binary <- function(stem, label, column) {
  if (!column %in% names(genes)) {
    add_missing(paste0(stem, "_violin/hist"), paste("missing columns:", column))
    return(NULL)
  }

  df <- genes %>%
    transmute(sample, value = .metric, present = coerce_binary_plot(.data[[column]])) %>%
    filter(!is.na(present), is.finite(value)) %>%
    mutate(present = factor(if_else(present, "Present", "Absent"), levels = c("Absent", "Present")))

  if (n_distinct(df$present) < 2L) {
    add_missing(paste0(stem, "_violin/hist"), "only one presence group")
    return(NULL)
  }

  legend_labels <- make_binary_labels(df)

  per_genome <- df %>%
    group_by(sample, present) %>%
    summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")
  p <- safe_wilcox_p(per_genome$median_value, per_genome$present)

  df_violin <- sample_per_group(df, "present", max_violin_rows_per_group)
  subtitle_violin <- wrap_text(
    paste0(
      "Violin drawn from a stratified sample for speed; p-value from per-genome medians. ",
      "All counts in labels are total genes before sampling."
    ),
    width = 88
  )

  p_violin <- ggplot(df_violin, aes(x = present, y = value, fill = present)) +
    geom_violin(trim = FALSE, alpha = 0.70, na.rm = TRUE, width = 0.95) +
    geom_boxplot(width = 0.10, outlier.alpha = 0.12, alpha = 0.85, na.rm = TRUE) +
    annotate_top_right(sig_label(p)) +
    scale_x_discrete(labels = legend_labels) +
    scale_fill_discrete(labels = legend_labels) +
    labs(
      x = NULL, y = metric, fill = NULL,
      title = paste(metric, "by", label),
      subtitle = subtitle_violin
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")

  p_hist <- ggplot(df, aes(x = value, fill = present)) +
    geom_histogram(aes(y = after_stat(count / 1e6)), bins = 60, alpha = 0.45, position = "identity", na.rm = TRUE) +
    annotate_top_right(sig_label(p)) +
    scale_fill_discrete(labels = legend_labels) +
    labs(
      x = metric, y = "Gene count (million)", fill = NULL,
      title = paste(metric, "distribution by", label),
      subtitle = wrap_text("Overlaid histograms use identical bins on the same x-axis. Counts in legend are total genes.", 90)
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")

  save_plot_pair(p_violin, paste0(stem, "_violin"), output_dir, 7.4, 5.5, formats)
  save_plot_pair(p_hist, paste0(stem, "_hist"), output_dir, 7.8, 5.2, formats)
  p_violin + labs(title = label, subtitle = subtitle_violin)
}

plot_continuous <- function(stem, label, column) {
  if (!column %in% names(genes)) {
    add_missing(stem, paste("missing columns:", column)); return(NULL)
  }
  df <- genes %>%
    transmute(sample, value = .metric, x = suppressWarnings(as.numeric(.data[[column]]))) %>%
    filter(is.finite(value), is.finite(x))
  if (nrow(df) < 100L || length(unique(df$x)) < 3L) {
    add_missing(stem, "too few finite values or nearly constant x"); return(NULL)
  }
  ann <- lm_annotation(df$x, df$value)
  if (nrow(df) > max_scatter_rows) {
    df_plot <- df %>% slice_sample(n = max_scatter_rows)
    subtitle <- paste0("Regression annotation is computed on all finite rows; ", max_scatter_rows, " points are sampled for plotting speed.")
  } else {
    df_plot <- df
    subtitle <- "Regression annotation is computed on all finite rows."
  }
  p <- ggplot(df_plot, aes(x = x, y = value)) +
    geom_point(alpha = 0.08, size = 0.45) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.8) +
    annotate_top_right(ann) +
    labs(
      x = label, y = metric,
      title = paste(metric, "vs", label),
      subtitle = wrap_text(subtitle, width = 90)
    ) +
    theme_minimal(base_size = 11)
  save_plot_pair(p, stem, output_dir, 7.0, 5.2, formats)
  p + labs(title = label, subtitle = wrap_text(subtitle, 60))
}

violin_panels <- purrr::pmap(binary_features, plot_binary)
violin_panels <- violin_panels[!vapply(violin_panels, is.null, logical(1))]

if (length(violin_panels) > 0L && requireNamespace("patchwork", quietly = TRUE)) {
  combined <- patchwork::wrap_plots(violin_panels, ncol = 2) +
    patchwork::plot_annotation(
      title = paste(metric, "by binary protein features"),
      subtitle = wrap_text("Each violin uses a bounded stratified sample for rendering only. Tests use per-genome medians; category labels show total gene counts before sampling.", 115)
    )
  save_plot_pair(combined, "binary_feature_violins_grid", output_dir, 12.0, 8.7, formats)
}

purrr::pwalk(continuous_features, plot_continuous)

if (nrow(missing_plots) == 0L) {
  missing_plots <- tibble(plot = "none", reason = "all requested plots generated")
}
readr::write_tsv(missing_plots, file.path(output_dir, "XXXXX_missing_plots.tsv"))
message("Done: ", output_dir)
