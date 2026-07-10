#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# plot_go_term_feature_suggestion_set.R
# -----------------------------------------------------------------------------
# Purpose
#   Recreate the GO-term subset plots from scripts_suggestions/usefull_scripts/
#   goterms_plots.R, using the canonical gene_features.tsv table.
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

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

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
    go_col <- intersect(c("go_id", "go_term", "term", "go_terms", "GO", "id"), names(tab))[[1]]
    if (!is.na(go_col)) return(extract_go_ids(tab[[go_col]]))
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
    add_missing(stem, paste("missing column:", column))
    return(invisible(NULL))
  }
  df <- subset_genes |>
    transmute(presence = coerce_binary(.data[[column]]), value = .metric) |>
    filter(!is.na(presence), is.finite(value)) |>
    mutate(presence = factor(if_else(presence, "Present", "Absent"), levels = c("Absent", "Present")))
  if (n_distinct(df$presence) < 2L) {
    add_missing(stem, "fewer than two presence groups")
    return(invisible(NULL))
  }
  p_box <- ggplot(df, aes(presence, value, fill = presence)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.35) +
    labs(x = label, y = metric, title = paste("Selected GO genes:", metric, "by", label)) +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  p_violin <- ggplot(df, aes(presence, value, fill = presence)) +
    geom_violin(trim = FALSE, alpha = 0.75) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.85) +
    labs(x = label, y = metric, title = paste("Selected GO genes:", metric, "by", label)) +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  save_plot_multi(p_box, paste0(stem, "_boxplot"), width = 7.2, height = 5.2)
  save_plot_multi(p_violin, paste0(stem, "_violin"), width = 7.2, height = 5.2)
}

plot_numeric_pair <- function(column, label, stem, trimmed = FALSE) {
  if (!column %in% names(subset_genes)) {
    add_missing(stem, paste("missing column:", column))
    return(invisible(NULL))
  }
  df <- subset_genes |>
    transmute(x = suppressWarnings(as.numeric(.data[[column]])), value = .metric) |>
    filter(is.finite(x), is.finite(value))
  if (nrow(df) < 20L || length(unique(df$x)) < 2L) {
    add_missing(stem, "too few finite non-constant values")
    return(invisible(NULL))
  }
  if (trimmed) {
    hi <- quantile(df$x, 0.99, na.rm = TRUE)
    df <- df |> filter(x <= hi)
  }
  p <- ggplot(df, aes(x, value)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_smooth(method = "loess", se = FALSE, formula = y ~ x) +
    labs(x = label, y = metric, title = paste("Selected GO genes:", metric, "vs", label)) +
    theme_minimal(base_size = 11)
  save_plot_multi(p, stem, width = 7.2, height = 5.2)
}

# Old plot names preserved where possible.
plot_binary_pair("pfam_present", "PFAM/domain present", "GOterms_domain")
plot_binary_pair("lcr_present", "LCR present", "GOterms_LCRpres")
plot_binary_pair("signal_peptide_present", "Signal peptide present", "GOterms_signal_presence")
plot_binary_pair("tm_present", "TM present", "GOterms_TM_presence")
plot_binary_pair("pfam_lcr_present", "PFAM-LCR overlap", "GOterms_PFAM_LCR_presence") # XXXXX if absent.
plot_numeric_pair("lcr_total_length", "Total LCR length [aa]", "GOterms_LCRlen_scatter")
plot_numeric_pair("lcr_total_length", "Total LCR length [aa]", "GOterms_LCRlen_TRIMMED_scatter", trimmed = TRUE)
plot_numeric_pair("tm_total_length", "Total TM length [aa]", "GOterms_TM_len_scatter")

readr::write_tsv(tibble(go_id = chosen_go), file.path(output_dir, "selected_go_terms_used.tsv"))
if (nrow(missing_plots) == 0L) missing_plots <- tibble(plot = "none", reason = "all possible GO-term feature plots generated")
readr::write_tsv(missing_plots, file.path(output_dir, "XXXXX_missing_plots.tsv"))
message("Done. GO-term subset plots written to: ", output_dir)
