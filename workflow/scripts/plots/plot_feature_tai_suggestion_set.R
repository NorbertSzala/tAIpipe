#!/usr/bin/env Rscript
# Fast, report-oriented feature-vs-tAI plots.
#
# Design choices:
#   - The binary-feature overview is owned by the canonical gene-feature script;
#     this legacy layer keeps only the non-redundant continuous-feature panels.
#   - Violin plots use bounded stratified sampling for rendering only; tests are
#     computed from per-genome medians, not from the sampled gene rows.
#   - Per-genome summaries are the unit of inference to avoid significance
#     inflation caused by millions of non-independent genes.
#   - Feature-specific sample sizes are shown in facet labels because one global
#     legend cannot represent different counts for every feature honestly.

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
# Individual binary violins/histograms duplicated the combined figure and are
# intentionally removed. The faceted split-violin grid is the canonical output.
legacy_binary_stems <- c("domains", "LCR_presence", "signal_presence", "TM_presence", "PFAM_LCR")
unlink(unlist(lapply(legacy_binary_stems, function(z) file.path(output_dir, paste0(z, c("_hist.png", "_hist.pdf", "_violin.png", "_violin.pdf"))))), force = TRUE)
unlink(file.path(output_dir, c("binary_feature_violins_grid.png", "binary_feature_violins_grid.pdf")), force = TRUE)

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

collect_binary <- function(stem, label, column) {
  if (!column %in% names(genes)) {
    add_missing(paste0(stem, "_split_violin"), paste("missing column:", column))
    return(tibble())
  }
  df <- genes %>%
    transmute(sample, value = .metric, present = coerce_binary_plot(.data[[column]])) %>%
    filter(!is.na(present), is.finite(value)) %>%
    mutate(present = factor(if_else(present, "Present", "Absent"), levels = c("Absent", "Present")))
  if (n_distinct(df$present) < 2L) {
    add_missing(paste0(stem, "_split_violin"), "only one presence group")
    return(tibble())
  }
  counts <- df %>% count(present, name = "n")
  absent_n <- counts$n[counts$present == "Absent"] %||% 0L
  present_n <- counts$n[counts$present == "Present"] %||% 0L
  per_genome <- df %>% group_by(sample, present) %>%
    summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")
  pval <- safe_paired_wilcox_p(
    per_genome, sample, present, median_value,
    level_order = c("Absent", "Present")
  )
  sample_per_group(df, "present", max_violin_rows_per_group) %>%
    mutate(
      feature = label,
      facet_label = paste0(label, "\nAbsent (", format(absent_n, big.mark = ","),
                           "); Present (", format(present_n, big.mark = ","), ")",
                           "\n", sig_label(pval))
    )
}

message(
  "Binary feature overview omitted here: the canonical gene_features/",
  "tai_by_binary_features plot contains the same comparison."
)

plot_continuous <- function(stem, label, column, max_x = NULL, suffix = NULL) {
  if (!column %in% names(genes)) {
    add_missing(stem, paste("missing columns:", column)); return(NULL)
  }
  df <- genes %>%
    transmute(sample, value = .metric, x = suppressWarnings(as.numeric(.data[[column]]))) %>%
    filter(is.finite(value), is.finite(x)) %>%
    filter(if (column %in% c("lcr_total_length", "tm_total_length", "lcr_count", "tm_count")) x > 0 else TRUE)
  if (!is.null(max_x)) df <- df %>% filter(x <= max_x)
  if (nrow(df) < 100L || length(unique(df$x)) < 3L) {
    add_missing(paste0(stem, suffix %||% ""), "too few finite values or nearly constant x"); return(NULL)
  }
  ann <- lm_annotation(df$x, df$value)
  if (nrow(df) > max_scatter_rows) {
    df_plot <- df %>% slice_sample(n = max_scatter_rows)
  } else {
    df_plot <- df
  }
  title_suffix <- if (is.null(max_x)) "" else paste0(" (x <= ", max_x, ")")
  p <- ggplot(df_plot, aes(x = x, y = value)) +
    geom_point(alpha = 0.08, size = 0.45, colour = "#1F4E79") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.8, colour = "grey25") +
    annotate_top_right(ann, size = 3.9) +
    labs(
      x = label,
      y = metric_axis_label(metric),
      title = paste0(label, title_suffix),
      subtitle = NULL
    ) +
    theme_minimal(base_size = 12.5) +
    theme(
      axis.text = element_text(size = 11.5),
      axis.title = element_text(size = 12.5),
      plot.title = element_text(size = 15.5, face = "bold")
    )
  p_full <- add_single_marginal_densities(
    p, df_plot, x, value, colour = "#1F4E79",
    top_height = 1.3, right_width = 1.3
  )
  out_stem <- paste0(stem, suffix %||% "")
  save_plot_pair(p_full, out_stem, output_dir, 9.0, 7.0, formats)
  invisible(p_full)
}


purrr::pwalk(continuous_features, plot_continuous)
plot_continuous("LCR_length", "Total LCR length [aa]", "lcr_total_length", max_x = 500, suffix = "_max500")
plot_continuous("LCR_number", "LCR count", "lcr_count", max_x = 20, suffix = "_max20")
plot_continuous("protein_length", "Protein length [aa]", "protein_length_aa", max_x = 3000, suffix = "_max3000")
plot_continuous("TM_length", "Total TM length [aa]", "tm_total_length", max_x = 300, suffix = "_max300")

if (nrow(missing_plots) == 0L) {
  missing_plots <- tibble(plot = "none", reason = "all requested plots generated")
}
readr::write_tsv(missing_plots, file.path(output_dir, "XXXXX_missing_plots.tsv"))
message("Done: ", output_dir)
