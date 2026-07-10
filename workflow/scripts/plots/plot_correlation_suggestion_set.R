#!/usr/bin/env Rscript
# Correlation heatmaps for the script_suggestions-compatible plot layer.
# The plotted variables are taken from script_suggestion_plots.*_correlation_variables
# in the active Snakemake configuration. Missing/non-informative variables are
# reported explicitly instead of being silently replaced by hard-coded defaults.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})
source("workflow/scripts/lib/plot_style_helpers.R")

parse_vector <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return(character())
  }
  values <- unlist(strsplit(paste(as.character(x), collapse = ","), ",", fixed = TRUE), use.names = FALSE)
  values <- str_trim(values)
  unique(values[nzchar(values)])
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    gene_table = "results/tables/gene_features.tsv",
    genome_table = "results/tables/genome_summary.tsv",
    output_dir = "results/plots/script_suggestions/correlations",
    method = "spearman",
    formats = "png,pdf",
    gene_variables = "tAI,CAI,ENC,GC,GC3s,protein_length_aa,lcr_total_length,tm_total_length",
    genome_variables = "mean_tAI,mean_CAI,mean_ENC,genome_gc,mean_GC3s,fraction_lcr,fraction_tm,fraction_signal_peptide,fraction_pfam"
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
    gene_table = snakemake@input[["gene_table"]] %||% snakemake@input[["gene_features"]],
    genome_table = snakemake@input[["genome_summary"]] %||% snakemake@input[["genome_table"]],
    output_dir = snakemake@params[["output_dir"]] %||% "results/plots/script_suggestions/correlations",
    method = snakemake@params[["method"]] %||% "spearman",
    formats = paste(snakemake@params[["formats"]] %||% c("png", "pdf"), collapse = ","),
    gene_variables = paste(snakemake@params[["gene_variables"]] %||% character(), collapse = ","),
    genome_variables = paste(snakemake@params[["genome_variables"]] %||% character(), collapse = ",")
  )
} else parse_args()

formats <- parse_vector(args$formats)
if (length(formats) == 0L) {
  formats <- c("png", "pdf")
}

method <- tolower(as.character(args$method))
if (!method %in% c("pearson", "spearman", "kendall")) {
  stop("Unsupported correlation method: ", args$method)
}

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

gene_vars_configured <- parse_vector(args$gene_variables)
genome_vars_configured <- parse_vector(args$genome_variables)

variable_status <- function(df, requested, level) {
  if (length(requested) == 0L) {
    return(tibble(level = character(), variable = character(), status = character(), n_finite = integer(), n_unique = integer()))
  }

  tibble(level = level, variable = requested) %>%
    mutate(
      present = variable %in% names(df),
      n_finite = vapply(variable, function(v) {
        if (!v %in% names(df)) return(0L)
        x <- suppressWarnings(as.numeric(df[[v]]))
        sum(is.finite(x))
      }, integer(1)),
      n_unique = vapply(variable, function(v) {
        if (!v %in% names(df)) return(0L)
        x <- suppressWarnings(as.numeric(df[[v]]))
        length(unique(x[is.finite(x)]))
      }, integer(1)),
      status = case_when(
        !present ~ "missing_column",
        n_finite < 3L ~ "too_few_finite_values",
        n_unique < 2L ~ "constant_or_single_value",
        TRUE ~ "used"
      )
    ) %>%
    select(level, variable, status, n_finite, n_unique)
}

cor_long <- function(df, vars, method) {
  status <- variable_status(df, vars, level = NA_character_)
  usable <- status %>% filter(status == "used") %>% pull(variable)
  if (length(usable) < 2L) {
    return(tibble())
  }

  mat <- df %>% transmute(across(all_of(usable), ~ suppressWarnings(as.numeric(.x))))
  cm <- suppressWarnings(cor(mat, use = "pairwise.complete.obs", method = method))

  as.data.frame(as.table(cm)) %>%
    as_tibble() %>%
    rename(variable_x = Var1, variable_y = Var2, correlation = Freq) %>%
    mutate(
      variable_x_label = clean_plot_label(variable_x),
      variable_y_label = clean_plot_label(variable_y)
    )
}

plot_cor <- function(d, title, stem, width = 8, height = 7) {
  if (nrow(d) == 0L) return(invisible(NULL))

  p <- ggplot(d, aes(x = variable_x_label, y = variable_y_label, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, fill = paste0(method, " r"), title = title) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11)
    )

  save_plot_pair(p, stem, args$output_dir, width, height, formats)
}

genes <- readr::read_tsv(args$gene_table, show_col_types = FALSE, progress = FALSE)
genomes <- readr::read_tsv(args$genome_table, show_col_types = FALSE, progress = FALSE)

gene_status <- variable_status(genes, gene_vars_configured, "gene")
genome_status <- variable_status(genomes, genome_vars_configured, "genome")
readr::write_tsv(bind_rows(gene_status, genome_status), file.path(args$output_dir, "correlation_variable_status.tsv"))

gene_vars <- gene_status %>% filter(status == "used") %>% pull(variable)
genome_vars <- genome_status %>% filter(status == "used") %>% pull(variable)

gene_cor <- cor_long(genes, gene_vars, method)
genome_cor <- cor_long(genomes, genome_vars, method)

plot_cor(gene_cor, "Gene-level correlations", "gene_level_correlation_heatmap", 9, 8)
plot_cor(genome_cor, "Genome-level correlations", "genome_level_correlation_heatmap", 9, 8)

readr::write_tsv(gene_cor, file.path(args$output_dir, "gene_level_correlation_heatmap_data.tsv"))
readr::write_tsv(genome_cor, file.path(args$output_dir, "genome_level_correlation_heatmap_data.tsv"))

message("Saved script_suggestions correlation outputs to: ", args$output_dir)
