#!/usr/bin/env Rscript
# Compute conservative gene-feature effects using genome as the unit of inference.
#
# This script replaces/augments pooled gene-level mixed models when very large
# gene counts make tiny p-values biologically uninformative. For each binary
# feature it estimates an adjusted within-genome effect, then tests whether
# those effects are consistently different from zero across genomes.
#
# Example:
# Rscript workflow/scripts/statistics/compute_gene_feature_effects_per_genome.R \
#   --gene-table results/tables/gene_features.tsv \
#   --features signal_peptide_present,tm_present,lcr_present,pfam_present \
#   --covariates log_protein_length_aa,GC3s \
#   --per-genome-output results/statistics/gene_feature_per_genome_effects.tsv \
#   --meta-output results/statistics/gene_feature_effect_meta_tests.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    key <- sub("^--", "", key)
    if (i == length(argv) || startsWith(argv[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- argv[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

get_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (is.null(value)) {
    if (required) stop("Missing required argument --", name)
    return(default)
  }
  value
}

as_integer_arg <- function(args, name, default) {
  value <- as.integer(get_arg(args, name, default))
  if (!is.finite(value)) stop("Argument --", name, " must be integer")
  value
}

args_raw <- parse_args()
args <- list(
  gene_table = get_arg(args_raw, "gene-table", required = TRUE),
  features = get_arg(args_raw, "features", required = TRUE),
  covariates = get_arg(args_raw, "covariates", "log_protein_length_aa,GC3s"),
  value_column = get_arg(args_raw, "value-column", "tAI"),
  sample_column = get_arg(args_raw, "sample-column", "sample"),
  min_genes = as_integer_arg(args_raw, "min-genes", 100L),
  min_class_genes = as_integer_arg(args_raw, "min-class-genes", 10L),
  min_genomes = as_integer_arg(args_raw, "min-genomes", 5L),
  fdr_method = get_arg(args_raw, "fdr-method", "BH"),
  per_genome_output = get_arg(args_raw, "per-genome-output", required = TRUE),
  meta_output = get_arg(args_raw, "meta-output", required = TRUE)
)

split_argument <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(character())
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  unique(values[nzchar(values)])
}

as_binary <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) {
    finite_values <- unique(x[is.finite(x)])
    if (all(finite_values %in% c(0, 1))) return(as.integer(x))
  }
  normalized <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(normalized))
  out[normalized %in% c("1", "true", "t", "yes", "y", "present")] <- 1L
  out[normalized %in% c("0", "false", "f", "no", "n", "absent")] <- 0L
  out
}

safe_zscore <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  finite <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(finite) < 2L) return(out)
  s <- stats::sd(x[finite])
  if (!is.finite(s) || s == 0) return(out)
  out[finite] <- (x[finite] - mean(x[finite])) / s
  out
}

numeric_has_variation <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  length(x) >= 2L && length(unique(x)) >= 2L
}

empty_per_genome <- function(feature, sample_id, status, n_genes = 0L,
                             n_present = 0L, n_absent = 0L,
                             used_covariates = character(), formula = NA_character_) {
  tibble(
    feature = feature,
    sample = sample_id,
    estimate = NA_real_,
    std_error = NA_real_,
    statistic = NA_real_,
    p_value_lm = NA_real_,
    mean_present = NA_real_,
    mean_absent = NA_real_,
    mean_delta_unadjusted = NA_real_,
    median_present = NA_real_,
    median_absent = NA_real_,
    median_delta_unadjusted = NA_real_,
    n_genes = as.integer(n_genes),
    n_present = as.integer(n_present),
    n_absent = as.integer(n_absent),
    used_covariates = paste(used_covariates, collapse = ";"),
    model_formula = formula,
    status = status
  )
}

fit_one_sample_feature <- function(sample_data, feature, requested_covariates,
                                   sample_column, value_column,
                                   min_genes, min_class_genes) {
  sample_id <- unique(sample_data[[sample_column]])[[1]]

  if (!feature %in% names(sample_data)) {
    return(empty_per_genome(feature, sample_id, "missing_feature_column"))
  }

  available_covariates <- intersect(requested_covariates, names(sample_data))
  if ("log_protein_length_aa" %in% requested_covariates &&
      !"log_protein_length_aa" %in% names(sample_data) &&
      "protein_length_aa" %in% names(sample_data)) {
    sample_data$log_protein_length_aa <- log1p(suppressWarnings(as.numeric(sample_data$protein_length_aa)))
    available_covariates <- union(available_covariates, "log_protein_length_aa")
  }

  model_data <- sample_data |>
    transmute(
      tAI_z = safe_zscore(.data[[value_column]]),
      feature_value = as_binary(.data[[feature]]),
      across(all_of(available_covariates), ~ suppressWarnings(as.numeric(.x)))
    ) |>
    filter(is.finite(tAI_z), !is.na(feature_value))

  n_present <- sum(model_data$feature_value == 1L, na.rm = TRUE)
  n_absent <- sum(model_data$feature_value == 0L, na.rm = TRUE)

  if (nrow(model_data) < min_genes) {
    return(empty_per_genome(feature, sample_id, "too_few_genes", nrow(model_data), n_present, n_absent))
  }
  if (n_present < min_class_genes || n_absent < min_class_genes) {
    return(empty_per_genome(feature, sample_id, "too_few_genes_in_one_class", nrow(model_data), n_present, n_absent))
  }

  usable_covariates <- available_covariates[vapply(
    model_data[available_covariates], numeric_has_variation, logical(1)
  )]

  keep_cols <- c("tAI_z", "feature_value", usable_covariates)
  model_data <- model_data |>
    select(all_of(keep_cols)) |>
    filter(if_all(everything(), ~ !is.na(.x)))

  if (length(usable_covariates) > 0L) {
    model_data <- model_data |>
      filter(if_all(all_of(usable_covariates), is.finite))
  }

  n_present <- sum(model_data$feature_value == 1L, na.rm = TRUE)
  n_absent <- sum(model_data$feature_value == 0L, na.rm = TRUE)

  if (nrow(model_data) < min_genes || n_present < min_class_genes || n_absent < min_class_genes) {
    return(empty_per_genome(
      feature, sample_id, "too_few_complete_cases", nrow(model_data), n_present, n_absent, usable_covariates
    ))
  }

  formula_text <- paste(
    "tAI_z ~ feature_value",
    if (length(usable_covariates) > 0L) paste("+", paste(usable_covariates, collapse = " + ")) else ""
  )

  fit <- tryCatch(
    stats::lm(stats::as.formula(formula_text), data = model_data),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(empty_per_genome(
      feature, sample_id, paste0("lm_error: ", conditionMessage(fit)),
      nrow(model_data), n_present, n_absent, usable_covariates, formula_text
    ))
  }

  coef_table <- summary(fit)$coefficients
  if (!"feature_value" %in% rownames(coef_table)) {
    return(empty_per_genome(
      feature, sample_id, "feature_term_not_estimable",
      nrow(model_data), n_present, n_absent, usable_covariates, formula_text
    ))
  }

  present_values <- model_data$tAI_z[model_data$feature_value == 1L]
  absent_values <- model_data$tAI_z[model_data$feature_value == 0L]

  tibble(
    feature = feature,
    sample = sample_id,
    estimate = unname(coef_table["feature_value", "Estimate"]),
    std_error = unname(coef_table["feature_value", "Std. Error"]),
    statistic = unname(coef_table["feature_value", "t value"]),
    p_value_lm = unname(coef_table["feature_value", "Pr(>|t|)"]),
    mean_present = mean(present_values, na.rm = TRUE),
    mean_absent = mean(absent_values, na.rm = TRUE),
    mean_delta_unadjusted = mean(present_values, na.rm = TRUE) - mean(absent_values, na.rm = TRUE),
    median_present = stats::median(present_values, na.rm = TRUE),
    median_absent = stats::median(absent_values, na.rm = TRUE),
    median_delta_unadjusted = stats::median(present_values, na.rm = TRUE) - stats::median(absent_values, na.rm = TRUE),
    n_genes = as.integer(nrow(model_data)),
    n_present = as.integer(n_present),
    n_absent = as.integer(n_absent),
    used_covariates = paste(usable_covariates, collapse = ";"),
    model_formula = formula_text,
    status = "ok"
  )
}

signed_rank_p <- function(x, min_n) {
  x <- x[is.finite(x)]
  if (length(x) < min_n) return(NA_real_)
  if (length(unique(x)) < 2L) return(NA_real_)
  tryCatch(
    stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

t_test_p <- function(x, min_n) {
  x <- x[is.finite(x)]
  if (length(x) < min_n) return(NA_real_)
  if (stats::sd(x) == 0) return(NA_real_)
  tryCatch(stats::t.test(x, mu = 0)$p.value, error = function(e) NA_real_)
}

summarize_feature <- function(per_genome, feature, min_genomes) {
  x <- per_genome |>
    filter(.data$feature == feature, .data$status == "ok", is.finite(.data$estimate)) |>
    pull(.data$estimate)

  if (length(x) < min_genomes) {
    return(tibble(
      feature = feature,
      analysis = "per_genome_effect_meta_test",
      n_genomes = length(x),
      mean_effect = if (length(x) > 0L) mean(x) else NA_real_,
      median_effect = if (length(x) > 0L) stats::median(x) else NA_real_,
      sd_effect = if (length(x) > 1L) stats::sd(x) else NA_real_,
      q25_effect = if (length(x) > 0L) unname(stats::quantile(x, 0.25, na.rm = TRUE)) else NA_real_,
      q75_effect = if (length(x) > 0L) unname(stats::quantile(x, 0.75, na.rm = TRUE)) else NA_real_,
      fraction_positive = if (length(x) > 0L) mean(x > 0) else NA_real_,
      fraction_negative = if (length(x) > 0L) mean(x < 0) else NA_real_,
      p_value_signed_rank = NA_real_,
      p_value_t_test = NA_real_,
      status = "too_few_genomes"
    ))
  }

  tibble(
    feature = feature,
    analysis = "per_genome_effect_meta_test",
    n_genomes = length(x),
    mean_effect = mean(x),
    median_effect = stats::median(x),
    sd_effect = stats::sd(x),
    q25_effect = unname(stats::quantile(x, 0.25, na.rm = TRUE)),
    q75_effect = unname(stats::quantile(x, 0.75, na.rm = TRUE)),
    fraction_positive = mean(x > 0),
    fraction_negative = mean(x < 0),
    p_value_signed_rank = signed_rank_p(x, min_genomes),
    p_value_t_test = t_test_p(x, min_genomes),
    status = "ok"
  )
}

features <- split_argument(args$features)
covariates <- split_argument(args$covariates)

if (length(features) == 0L) stop("No features supplied")
if (!file.exists(args$gene_table)) stop("Missing --gene-table: ", args$gene_table)

genes <- read_tsv(args$gene_table, show_col_types = FALSE, progress = FALSE)
required <- c(args$sample_column, args$value_column)
missing_required <- setdiff(required, names(genes))
if (length(missing_required) > 0L) {
  stop("Gene table lacks required columns: ", paste(missing_required, collapse = ", "))
}

missing_features <- setdiff(features, names(genes))
if (length(missing_features) > 0L) {
  warning("Missing feature columns will be reported as missing: ", paste(missing_features, collapse = ", "))
}

per_genome <- map_dfr(features, function(feature) {
  genes |>
    filter(!is.na(.data[[args$sample_column]]), nzchar(trimws(as.character(.data[[args$sample_column]])))) |>
    group_split(.data[[args$sample_column]]) |>
    map_dfr(~ fit_one_sample_feature(
      sample_data = .x,
      feature = feature,
      requested_covariates = covariates,
      sample_column = args$sample_column,
      value_column = args$value_column,
      min_genes = args$min_genes,
      min_class_genes = args$min_class_genes
    ))
})

meta <- map_dfr(features, ~ summarize_feature(per_genome, .x, args$min_genomes)) |>
  mutate(
    q_value_signed_rank = p.adjust(p_value_signed_rank, method = args$fdr_method),
    q_value_t_test = p.adjust(p_value_t_test, method = args$fdr_method)
  )

dir.create(dirname(args$per_genome_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$meta_output), recursive = TRUE, showWarnings = FALSE)
write_tsv(per_genome, args$per_genome_output, na = "NA")
write_tsv(meta, args$meta_output, na = "NA")

message("Saved per-genome effects: ", args$per_genome_output)
message("Saved meta-tests: ", args$meta_output)
