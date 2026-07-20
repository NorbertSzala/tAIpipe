#!/usr/bin/env Rscript
# Fast statistical summary for tAIpipe.
#
# This replaces the legacy all-gene mixed model with genome-level summaries for
# binary feature effects. The biological replication unit is the genome, not the
# individual gene, which avoids extremely slow lme4 fits on millions of genes and
# reduces pseudo-replication in gene-level feature tests.
#
# Outputs keep the old column names used by downstream plotting scripts where
# possible. The analysis label is changed to describe the new method.

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

parser <- ArgumentParser(description = "Compute tAIpipe statistical tests using genome-level summaries")
parser$add_argument("--gene-table", required = TRUE)
parser$add_argument("--genome-table", required = TRUE)
parser$add_argument("--gene-feature-output", required = TRUE)
parser$add_argument("--genome-group-output", required = TRUE)
parser$add_argument("--binary-features", default = "signal_peptide_present,tm_present,lcr_present,pfam_present")
parser$add_argument("--gene-covariates", default = "log_protein_length_aa,GC3s")
parser$add_argument("--genome-metrics", default = "mean_tAI,mean_CAI,mean_GC,mean_GC3s,mean_ENC,mean_delta_ENC")
parser$add_argument("--group-variables", default = "phylum,lifestyle")
parser$add_argument("--fdr-method", default = "BH")
parser$add_argument("--min-genes-per-sample", type = "integer", default = 100L)
parser$add_argument("--min-class-genes", type = "integer", default = 10L)
parser$add_argument("--min-genomes", type = "integer", default = 5L)
args <- parser$parse_args()

split_argument <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(character())
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  unique(values[nzchar(values)])
}

safe_zscore <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  finite <- is.finite(x)
  result <- rep(NA_real_, length(x))
  if (sum(finite) < 2L) return(result)
  spread <- stats::sd(x[finite])
  if (!is.finite(spread) || spread == 0) return(result)
  result[finite] <- (x[finite] - mean(x[finite])) / spread
  result
}

as_binary <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) {
    valid <- unique(x[is.finite(x)])
    if (all(valid %in% c(0, 1))) return(as.integer(x))
  }
  normalized <- tolower(trimws(as.character(x)))
  result <- rep(NA_integer_, length(normalized))
  result[normalized %in% c("1", "true", "t", "yes", "y", "present")] <- 1L
  result[normalized %in% c("0", "false", "f", "no", "n", "absent")] <- 0L
  result
}

empty_gene_result <- function(feature, status, n_genes = 0L, n_genomes = 0L) {
  tibble(
    analysis = "gene_feature_per_genome_median_difference",
    feature = feature,
    term = "present_minus_absent",
    estimate = NA_real_,
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    n_genes = as.integer(n_genes),
    n_genomes = as.integer(n_genomes),
    covariates = "not_used_in_fast_per_genome_summary",
    model_formula = "per-genome median(tAI_z | present) - median(tAI_z | absent)",
    status = status
  )
}

fit_binary_feature_per_genome <- function(gene_data, feature_name) {
  if (!feature_name %in% names(gene_data)) {
    return(empty_gene_result(feature_name, "missing_feature_column"))
  }

  d <- gene_data %>%
    transmute(
      sample = as.character(sample),
      tAI_z = suppressWarnings(as.numeric(tAI_z)),
      feature_value = as_binary(.data[[feature_name]])
    ) %>%
    filter(!is.na(sample), nzchar(sample), is.finite(tAI_z), !is.na(feature_value))

  if (nrow(d) < 20L || dplyr::n_distinct(d$sample) < args$min_genomes ||
      dplyr::n_distinct(d$feature_value) != 2L) {
    return(empty_gene_result(
      feature_name, "insufficient_data", nrow(d), dplyr::n_distinct(d$sample)
    ))
  }

  per_genome <- d %>%
    group_by(sample, feature_value) %>%
    summarise(
      n = dplyr::n(),
      median_tAI_z = median(tAI_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n >= args$min_class_genes) %>%
    mutate(feature_value = if_else(feature_value == 1L, "present", "absent")) %>%
    tidyr::pivot_wider(
      id_cols = sample,
      names_from = feature_value,
      values_from = c(n, median_tAI_z)
    ) %>%
    filter(
      !is.na(median_tAI_z_present),
      !is.na(median_tAI_z_absent),
      !is.na(n_present),
      !is.na(n_absent)
    ) %>%
    mutate(effect = median_tAI_z_present - median_tAI_z_absent)

  n_genomes <- nrow(per_genome)
  if (n_genomes < args$min_genomes) {
    return(empty_gene_result(feature_name, "insufficient_informative_genomes", nrow(d), n_genomes))
  }

  # Compute the signed-rank p-value independently of the CI. Requesting a
  # confidence interval can fail for degenerate/tie-heavy data even when the
  # signed-rank statistic itself is still available.
  test <- tryCatch(
    stats::wilcox.test(per_genome$effect, mu = 0, exact = FALSE),
    error = function(e) e
  )
  if (inherits(test, "error")) {
    return(empty_gene_result(feature_name, paste0("wilcox_error: ", conditionMessage(test)), nrow(d), n_genomes))
  }

  ci_test <- tryCatch(
    suppressWarnings(stats::wilcox.test(
      per_genome$effect, mu = 0, exact = FALSE,
      conf.int = TRUE, conf.level = 0.95
    )),
    error = function(e) NULL
  )

  # Wilcoxon with conf.int=TRUE returns the Hodges-Lehmann pseudomedian and its
  # 95% CI. The previous implementation stored the 2.5th/97.5th percentiles of
  # genome effects in conf_low/conf_high; those described the spread of genome
  # effects and were not confidence limits for the central effect.
  est <- if (!is.null(ci_test$estimate) && is.finite(unname(ci_test$estimate))) {
    unname(ci_test$estimate)
  } else {
    stats::median(per_genome$effect, na.rm = TRUE)
  }
  ci <- if (!is.null(ci_test$conf.int) && length(ci_test$conf.int) == 2L) {
    unname(ci_test$conf.int)
  } else {
    c(NA_real_, NA_real_)
  }
  se <- NA_real_  # no conventional standard error is reported for this non-parametric pseudomedian

  tibble(
    analysis = "gene_feature_per_genome_median_difference",
    feature = feature_name,
    term = "present_minus_absent",
    estimate = est,
    std_error = se,
    statistic = unname(test$statistic),
    p_value = test$p.value,
    conf_low = ci[[1]],
    conf_high = ci[[2]],
    n_genes = as.integer(nrow(d)),
    n_genomes = as.integer(n_genomes),
    covariates = "not_used_in_fast_per_genome_summary",
    model_formula = "per-genome median(tAI_z | present) - median(tAI_z | absent); Hodges-Lehmann pseudomedian with 95% Wilcoxon CI; signed-rank test against 0",
    status = "ok"
  )
}

empty_genome_result <- function(metric, group, status, n = 0L, k = 0L) {
  tibble(
    analysis = "genome_level_group_test",
    metric = metric,
    group_variable = group,
    test = NA_character_,
    comparison = NA_character_,
    effect = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    n_genomes = as.integer(n),
    n_groups = as.integer(k),
    status = status
  )
}

run_genome_group_test <- function(genome_data, metric_name, group_name) {
  if (!all(c(metric_name, group_name) %in% names(genome_data))) {
    return(empty_genome_result(metric_name, group_name, "missing_column"))
  }

  test_data <- genome_data %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[metric_name]])),
      group = trimws(as.character(.data[[group_name]]))
    ) %>%
    filter(is.finite(value), !is.na(group), nzchar(group)) %>%
    mutate(group = factor(group)) %>%
    droplevels()

  group_count <- nlevels(test_data$group)
  if (nrow(test_data) < 4L || group_count < 2L) {
    return(empty_genome_result(metric_name, group_name, "insufficient_data", nrow(test_data), group_count))
  }

  if (group_count == 2L) {
    if (any(table(test_data$group) < 2L)) {
      return(empty_genome_result(metric_name, group_name, "insufficient_group_size", nrow(test_data), group_count))
    }
    levels_group <- levels(test_data$group)
    result <- stats::wilcox.test(value ~ group, data = test_data, exact = FALSE)
    medians <- tapply(test_data$value, test_data$group, median)
    effect <- unname(medians[[2]] - medians[[1]])

    return(tibble(
      analysis = "genome_level_group_test",
      metric = metric_name,
      group_variable = group_name,
      test = "Wilcoxon rank-sum",
      comparison = paste(levels_group, collapse = " vs "),
      effect = effect,
      statistic = unname(result$statistic),
      p_value = result$p.value,
      n_genomes = nrow(test_data),
      n_groups = group_count,
      status = "ok"
    ))
  }

  result <- stats::kruskal.test(value ~ group, data = test_data)
  h_statistic <- unname(result$statistic)
  epsilon_squared <- max(0, (h_statistic - group_count + 1) / (nrow(test_data) - group_count))

  tibble(
    analysis = "genome_level_group_test",
    metric = metric_name,
    group_variable = group_name,
    test = "Kruskal-Wallis",
    comparison = "global",
    effect = epsilon_squared,
    statistic = h_statistic,
    p_value = result$p.value,
    n_genomes = nrow(test_data),
    n_groups = group_count,
    status = "ok"
  )
}

message("Reading gene table: ", args$gene_table)
gene_table <- readr::read_tsv(args$gene_table, show_col_types = FALSE, progress = FALSE)
message("Reading genome table: ", args$genome_table)
genome_table <- readr::read_tsv(args$genome_table, show_col_types = FALSE, progress = FALSE)

required_gene <- c("sample", "tAI")
missing_gene <- setdiff(required_gene, names(gene_table))
if (length(missing_gene) > 0L) stop("Gene table lacks columns: ", paste(missing_gene, collapse = ", "))

if (!"tAI_z" %in% names(gene_table)) {
  gene_table <- gene_table %>% group_by(sample) %>% mutate(tAI_z = safe_zscore(tAI)) %>% ungroup()
}

binary_features <- split_argument(args$binary_features)
if (length(binary_features) == 0L) stop("No binary features were configured")

message("Computing fast per-genome binary-feature effects")
gene_results <- purrr::map_dfr(binary_features, ~ fit_binary_feature_per_genome(gene_table, .x)) %>%
  mutate(q_value = p.adjust(p_value, method = args$fdr_method))

genome_metrics <- split_argument(args$genome_metrics)
group_variables <- split_argument(args$group_variables)
if (length(genome_metrics) == 0L || length(group_variables) == 0L) {
  stop("Genome metrics and group variables must not be empty")
}

message("Computing genome-level group tests")
genome_results <- tidyr::crossing(metric = genome_metrics, group_variable = group_variables) %>%
  purrr::pmap_dfr(~ run_genome_group_test(genome_table, ..1, ..2)) %>%
  mutate(q_value = p.adjust(p_value, method = args$fdr_method))

dir.create(dirname(args$gene_feature_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$genome_group_output), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(gene_results, args$gene_feature_output, na = "NA")
readr::write_tsv(genome_results, args$genome_group_output, na = "NA")
message("Saved: ", args$gene_feature_output)
message("Saved: ", args$genome_group_output)
