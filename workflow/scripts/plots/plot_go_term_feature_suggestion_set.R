#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# plot_go_term_feature_suggestion_set.R
# -----------------------------------------------------------------------------
# Purpose
#   Recreate the GO-term subset plots from scripts_suggestions/usefull_scripts/
#   goterms_plots.R as one combined feature grid using gene_features.tsv.
#
# Critical constraint
#   The original script uses a hand-made chosen_GOterms.tsv. The current project
#   does not define which GO terms should be selected for this figure. Therefore
#   this script does NOT guess. It generates GO-term feature plots only when an
#   explicit --chosen-go-terms file is provided.
#
# Accepted chosen-go-terms formats
#   1) One GO ID per line, e.g. GO:0006412
#   2) TSV/CSV with a column named go_id, go_term, term, or go_terms
#
# If chosen GO terms are missing or invalid, the script writes:
#   XXXXX_GO_TERMS_INPUT_REQUIRED.txt
# and exits without generating biological plots.
# -----------------------------------------------------------------------------

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
    chosen_go_terms = "XXXXX",
    output_dir = "results/plots/script_suggestions/go_term_features",
    metric = "tAI",
    formats = "png,pdf"
  )
  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--")) stop("Unexpected argument: ", args[[i]])
    if (i == length(args)) stop("Missing value for argument: ", args[[i]])
    key <- gsub("-", "_", sub("^--", "", args[[i]]))
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

if (exists("snakemake")) {
  args <- list(
    gene_table = snakemake@input[["gene_features"]],
    chosen_go_terms = snakemake@input[["chosen_go_terms"]] %||% "XXXXX",
    output_dir = snakemake@params[["output_dir"]] %||% "results/plots/script_suggestions/go_term_features",
    metric = snakemake@params[["metric"]] %||% "tAI",
    formats = paste(snakemake@params[["formats"]] %||% c("png", "pdf"), collapse = ",")
  )
} else {
  args <- parse_args()
}

output_dir <- args$output_dir
metric <- args$metric
formats <- str_split(args$formats, ",", simplify = FALSE)[[1]] |> trimws()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(output_dir, c(
  "GOterms_domain_boxplot.png", "GOterms_domain_boxplot.pdf",
  "GOterms_LCRpres_boxplot.png", "GOterms_LCRpres_boxplot.pdf",
  "GOterms_signal_presence_boxplot.png", "GOterms_signal_presence_boxplot.pdf",
  "GOterms_TM_presence_boxplot.png", "GOterms_TM_presence_boxplot.pdf",
  "GOterms_PFAM_LCR_presence_boxplot.png", "GOterms_PFAM_LCR_presence_boxplot.pdf",
  "GOterms_domain_violin.png", "GOterms_domain_violin.pdf",
  "GOterms_LCRpres_violin.png", "GOterms_LCRpres_violin.pdf",
  "GOterms_signal_presence_violin.png", "GOterms_signal_presence_violin.pdf",
  "GOterms_TM_presence_violin.png", "GOterms_TM_presence_violin.pdf",
  "GOterms_PFAM_LCR_presence_violin.png", "GOterms_PFAM_LCR_presence_violin.pdf",
  "GOterms_LCRlen_scatter.png", "GOterms_LCRlen_scatter.pdf",
  "GOterms_LCRlen_TRIMMED_scatter.png", "GOterms_LCRlen_TRIMMED_scatter.pdf",
  "GOterms_TM_len_scatter.png", "GOterms_TM_len_scatter.pdf"
)), force = TRUE)

write_xxxxx_and_exit <- function(reason) {
  path <- file.path(output_dir, "XXXXX_GO_TERMS_INPUT_REQUIRED.txt")
  writeLines(c(
    "XXXXX",
    "GO-term feature plots were not generated.",
    paste("Reason:", reason),
    "Provide an explicit --chosen-go-terms file with GO IDs to avoid arbitrary GO-term selection."
  ), path)
  message("XXXXX ", reason)
  quit(save = "no", status = 0)
}

if (!file.exists(args$gene_table)) stop("Missing gene table: ", args$gene_table)
if (identical(args$chosen_go_terms, "XXXXX") || !file.exists(args$chosen_go_terms)) {
  write_xxxxx_and_exit("missing explicit chosen GO-term list")
}

extract_go_ids <- function(x) {
  ids <- str_extract_all(as.character(x), regex("GO:[0-9]{7}", ignore_case = TRUE))
  ids <- unlist(ids, use.names = FALSE)
  unique(toupper(ids[!is.na(ids) & nzchar(ids)]))
}

read_chosen_go_terms <- function(path) {
  # Try delimited table first. If it fails or has no useful column, fall back to
  # reading raw lines.
  tab <- tryCatch(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(tab) && ncol(tab) > 0L) {
    candidates <- intersect(c("go_id", "go_term", "term", "go_terms", "GO", "id"), names(tab))
    if (length(candidates) > 0L) return(extract_go_ids(tab[[candidates[[1]]]]))
  }
  extract_go_ids(readLines(path, warn = FALSE))
}

chosen_go <- read_chosen_go_terms(args$chosen_go_terms)
if (length(chosen_go) == 0L) write_xxxxx_and_exit("chosen GO-term file contains no valid GO:0000000 identifiers")

genes <- readr::read_tsv(args$gene_table, show_col_types = FALSE, progress = FALSE)
required <- c("sample", "gene_id", "go_terms", metric)
missing <- setdiff(required, names(genes))
if (length(missing) > 0L) write_xxxxx_and_exit(paste("gene_features.tsv lacks columns:", paste(missing, collapse = ", ")))

# Filter genes annotated with any explicitly selected GO term.
go_pattern <- paste(chosen_go, collapse = "|")
subset_genes <- genes |>
  mutate(.metric = suppressWarnings(as.numeric(.data[[metric]]))) |>
  filter(is.finite(.metric), !is.na(go_terms), str_detect(go_terms, regex(go_pattern, ignore_case = TRUE)))

if (nrow(subset_genes) < 20L) {
  write_xxxxx_and_exit("fewer than 20 genes matched the selected GO terms")
}

save_plot_multi <- function(plot, stem, width = 8, height = 6) {
  for (fmt in formats) {
    ggsave(file.path(output_dir, paste0(stem, ".", fmt)), plot = plot, width = width, height = height, dpi = 300, limitsize = FALSE)
  }
}

coerce_binary <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x) || is.integer(x)) return(ifelse(is.na(x), NA, x > 0))
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "present", "presence")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "absent", "absence")] <- FALSE
  out
}

missing_plots <- tibble(plot = character(), reason = character())
add_missing <- function(plot, reason) {
  missing_plots <<- bind_rows(missing_plots, tibble(plot = plot, reason = reason))
  message("XXXXX skipped ", plot, ": ", reason)
}

plot_binary_pair <- function(column, label, stem) {
  if (!column %in% names(subset_genes)) {
    add_missing(stem, paste("missing column:", column)); return(invisible(NULL))
  }
  df <- subset_genes |>
    transmute(presence = coerce_binary(.data[[column]]), value = .metric) |>
    filter(!is.na(presence), is.finite(value)) |>
    mutate(presence = factor(if_else(presence, "Present", "Absent"), levels = c("Absent", "Present")))
  if (n_distinct(df$presence) < 2L) {
    add_missing(stem, "fewer than two presence groups"); return(invisible(NULL))
  }
  counts <- df |> count(presence, name = "n")
  counts <- counts %>%
    mutate(
      count_x = if_else(presence == "Absent", 0.82, 1.18),
      count_label = format(n, scientific = FALSE, trim = TRUE, big.mark = "")
    )
  df <- df %>%
    mutate(
      feature = label,
      split_side = split_side_from_level(presence, "Absent")
    )
  ggplot(
    df,
    aes(x = 1, y = value, fill = presence, group = presence, split_side = split_side)
  ) +
    geom_split_violin_project(trim = TRUE, alpha = 0.75, width = 0.94, scale = "width",
                              colour = "grey30", linewidth = 0.25) +
    geom_boxplot(
      aes(group = presence), width = 0.08, position = position_dodge(width = 0.18),
      outlier.shape = NA, alpha = 0.88, fill = "white"
    ) +
    geom_text(
      data = counts, aes(x = count_x, y = -Inf, label = count_label),
      inherit.aes = FALSE, vjust = 2.0, size = 2.7, colour = "grey25"
    ) +
    scale_x_continuous(breaks = NULL) +
    scale_fill_manual(values = binary_grey_values(levels(df$presence)), drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = metric_axis_label(metric), fill = NULL, title = label) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom", plot.title = element_text(size = 11.5, face = "bold", hjust = 0.5),
      plot.margin = margin(6, 6, 20, 6)
    )
}

plot_numeric_pair <- function(column, label, stem, trimmed = FALSE) {
  if (!column %in% names(subset_genes)) {
    add_missing(stem, paste("missing column:", column)); return(invisible(NULL))
  }
  df <- subset_genes |>
    transmute(x = suppressWarnings(as.numeric(.data[[column]])), value = .metric) |>
    filter(is.finite(x), x > 0, is.finite(value))
  if (nrow(df) < 20L || length(unique(df$x)) < 2L) {
    add_missing(stem, "too few positive, finite, non-constant values"); return(invisible(NULL))
  }
  if (trimmed) {
    hi <- quantile(df$x, 0.99, na.rm = TRUE)
    df <- df |> filter(x <= hi)
  }
  ann <- lm_annotation(df$x, df$value)
  p <- ggplot(df, aes(x, value)) +
    geom_point(alpha = 0.22, size = 0.55, colour = "#1F4E79") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey25", linewidth = 0.8) +
    annotate_top_right(ann, size = 3.1) +
    labs(x = label, y = metric_axis_label(metric), title = paste("Selected GO genes vs", label),
         subtitle = NULL) +
    theme_minimal(base_size = 12)
  p_full <- add_marginal_densities(
    p, df |> mutate(distribution = "Selected GO genes"), x, value, distribution,
    c("Selected GO genes" = "#1F4E79"), top_height = 1.35, right_width = 1.35
  )
  p_full
}

# Binary split violins and continuous relationships are collected into one grid
# so the legacy layer no longer emits seven partially redundant standalone files.
feature_plots <- Filter(Negate(is.null), list(
  plot_binary_pair("pfam_present", "PFAM/domain present", "GOterms_domain"),
  plot_binary_pair("lcr_present", "LCR present", "GOterms_LCRpres"),
  plot_binary_pair("signal_peptide_present", "Signal peptide present", "GOterms_signal_presence"),
  plot_binary_pair("tm_present", "TM present", "GOterms_TM_presence"),
  plot_binary_pair("pfam_lcr_present", "PFAM-LCR overlap", "GOterms_PFAM_LCR_presence"),
  plot_numeric_pair("lcr_total_length", "Total LCR length [aa]", "GOterms_LCRlen_scatter"),
  plot_numeric_pair("lcr_total_length", "Total LCR length [aa]", "GOterms_LCRlen_TRIMMED_scatter", trimmed = TRUE),
  plot_numeric_pair("tm_total_length", "Total TM length [aa]", "GOterms_TM_len_scatter")
))
if (length(feature_plots) > 0L) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for the combined GO-term feature grid.")
  }
  feature_grid <- patchwork::wrap_plots(feature_plots, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = "Selected GO genes: structural features and tAI",
      subtitle = "Split violins show absence (left, grey) and presence (right, blue); continuous panels show linear trends."
    ) &
    theme(legend.position = "bottom")
  save_plot_multi(feature_grid, "GOterms_feature_overview_grid", width = 15.5, height = 19.0)
} else {
  add_missing("GOterms_feature_overview_grid", "no feature panel could be generated")
}

readr::write_tsv(tibble(go_id = chosen_go), file.path(output_dir, "selected_go_terms_used.tsv"))
if (nrow(missing_plots) == 0L) missing_plots <- tibble(plot = "none", reason = "all possible GO-term feature plots generated")
readr::write_tsv(missing_plots, file.path(output_dir, "XXXXX_missing_plots.tsv"))
message("Done. GO-term subset plots written to: ", output_dir)
