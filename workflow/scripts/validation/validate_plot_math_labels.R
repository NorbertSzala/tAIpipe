#!/usr/bin/env Rscript

# Standalone validation of mathematical plot labels.
# This file must be executed with Rscript and does not require a Snakemake
# installation or the Snakemake-injected R object.

suppressPackageStartupMessages({
  library(stringr)
  library(latex2exp)
})

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) == 0L) {
    return(NA_character_)
  }

  normalizePath(
    sub("^--file=", "", file_arg[[1]]),
    winslash = "/",
    mustWork = TRUE
  )
}

script_path <- get_script_path()

if (is.na(script_path)) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  # workflow/scripts/validation/<this-file> -> repository root
  repo_root <- normalizePath(
    file.path(dirname(script_path), "..", "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
}

helper_path <- file.path(
  repo_root,
  "workflow",
  "scripts",
  "lib",
  "plot_style_helpers.R"
)

if (!file.exists(helper_path)) {
  stop(
    "Cannot locate plot_style_helpers.R at: ",
    helper_path,
    call. = FALSE
  )
}

helper_lines <- readLines(helper_path, warn = FALSE)

if (any(grepl("snakemake@", helper_lines, fixed = TRUE))) {
  stop(
    "plot_style_helpers.R unexpectedly contains a top-level Snakemake ",
    "reference. The wrong file may have been copied to: ",
    helper_path,
    call. = FALSE
  )
}

sys.source(helper_path, envir = .GlobalEnv)

required_functions <- c(
  "tex_label",
  "metric_axis_label",
  "metric_parsed_label",
  "metric_cv_label",
  "correlation_axis_labels"
)

missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]

if (length(missing_functions) > 0L) {
  stop(
    "Missing functions after sourcing plot_style_helpers.R: ",
    paste(missing_functions, collapse = ", "),
    call. = FALSE
  )
}

labels <- c(
  "$\\mathrm{tRNA}$ QC flags",
  "Binary features in $\\mathrm{tAI}$ distribution tails",
  "Continuous features in $\\mathrm{tAI}$ distribution tails",
  "$\\mathrm{GC}_{\\mathrm{genome}}$ and $\\overline{\\mathrm{GC}_{3s}}$ in Dikarya vs non-Dikarya",
  "Dashed line: $\\log_{2}(\\mathrm{OR})=0$; positive values indicate enrichment in the indicated $\\mathrm{tAI}$ tail",
  "Codon-level $\\mathrm{tRNA}$ profile per genome",
  "Cys-TGT-assigned $\\mathrm{tRNA}$s among all $\\mathrm{tRNA}$s [%]",
  "Selected GO terms from $\\mathrm{tAI}$-tail enrichment",
  "LCR length vs $\\mathrm{tAI}_{z}$ by PFAM-LCR overlap status; lowest 90% of lengths"
)

metrics <- c(
  "tAI", "tAI_z", "CAI", "ENC", "delta_ENC", "GC", "mean_GC",
  "genome_gc", "GC3s", "mean_GC3s", "genome_RSCU", "RSCU",
  "tRNA_weight", "trna_copy_number", "codon_frequency", "codon_count"
)

correlation_variables <- c(
  "tAI", "CAI", "ENC", "delta_ENC", "GC", "mean_GC", "GC3s",
  "protein_length_aa", "log_protein_length_aa", "lcr_total_length",
  "tm_total_length", "fraction_lcr", "fraction_tm"
)

for (label in labels) {
  tex_label(label)
  tex_label(label, output = "character")
}

for (metric in metrics) {
  metric_axis_label(metric)
  metric_parsed_label(metric)
  metric_cv_label(metric)
}

correlation_axis_labels(correlation_variables)

cat(
  "OK: validated ", length(labels), " mixed labels, ",
  length(metrics), " metric labels and ",
  length(correlation_variables), " correlation labels.\n",
  sep = ""
)
