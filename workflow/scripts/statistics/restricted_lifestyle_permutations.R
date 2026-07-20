#!/usr/bin/env Rscript
# Restricted permutation tests for genome-level metrics.
#
# Main use case in this project:
#   test whether lifestyle explains genome-level tAI/CAI/GC3s differences
#   after controlling for phylum, while permuting lifestyle labels only within phylum.
#
# Example:
# Rscript workflow/scripts/statistics/restricted_lifestyle_permutations.R \
#   --genome-table results/tables/genome_summary.tsv \
#   --responses mean_tAI,median_tAI,mean_CAI,median_CAI,mean_GC3s,median_GC3s \
#   --predictor lifestyle \
#   --strata phylum \
#   --n-perm 9999 \
#   --seed 1 \
#   --output results/statistics/lifestyle_within_phylum_permutation_tests.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
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

split_csv <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(character())
  x <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  unique(x[nzchar(x)])
}

as_numeric_arg <- function(args, name, default) {
  value <- as.numeric(get_arg(args, name, default))
  if (!is.finite(value)) stop("Argument --", name, " must be numeric")
  value
}

as_integer_arg <- function(args, name, default) {
  value <- as.integer(get_arg(args, name, default))
  if (!is.finite(value)) stop("Argument --", name, " must be integer")
  value
}

# F statistic for predictor adjusted for strata by comparing nested linear models.
# Full:    response ~ strata + predictor
# Reduced: response ~ strata
# This is less sensitive to formula order than extracting a row from anova(full).
observed_f_stat <- function(data, response, predictor, strata) {
  model_data <- data %>%
    transmute(
      response = suppressWarnings(as.numeric(.data[[response]])),
      predictor = as.factor(.data[[predictor]]),
      strata = as.factor(.data[[strata]])
    ) %>%
    filter(is.finite(response), !is.na(predictor), !is.na(strata)) %>%
    droplevels()

  if (nrow(model_data) < 5L) {
    return(list(f_stat = NA_real_, df_num = NA_real_, df_den = NA_real_, n = nrow(model_data), status = "too_few_rows"))
  }
  if (nlevels(model_data$predictor) < 2L) {
    return(list(f_stat = NA_real_, df_num = NA_real_, df_den = NA_real_, n = nrow(model_data), status = "predictor_has_one_level"))
  }
  if (nlevels(model_data$strata) < 2L) {
    return(list(f_stat = NA_real_, df_num = NA_real_, df_den = NA_real_, n = nrow(model_data), status = "strata_has_one_level"))
  }

  full_fit <- tryCatch(lm(response ~ strata + predictor, data = model_data), error = function(e) e)
  reduced_fit <- tryCatch(lm(response ~ strata, data = model_data), error = function(e) e)
  if (inherits(full_fit, "error") || inherits(reduced_fit, "error")) {
    return(list(f_stat = NA_real_, df_num = NA_real_, df_den = NA_real_, n = nrow(model_data), status = "lm_error"))
  }

  cmp <- tryCatch(anova(reduced_fit, full_fit), error = function(e) e)
  if (inherits(cmp, "error") || nrow(cmp) < 2L || !"F" %in% names(cmp)) {
    return(list(f_stat = NA_real_, df_num = NA_real_, df_den = NA_real_, n = nrow(model_data), status = "anova_error"))
  }

  list(
    f_stat = unname(cmp$F[[2L]]),
    df_num = unname(cmp$Df[[2L]]),
    df_den = unname(stats::df.residual(full_fit)),
    n = nrow(model_data),
    status = if (is.finite(cmp$F[[2L]])) "ok" else "nonfinite_f"
  )
}

permute_within_strata <- function(x, strata) {
  out <- x
  groups <- split(seq_along(x), strata, drop = TRUE)
  for (idx in groups) {
    if (length(idx) > 1L) {
      out[idx] <- sample(out[idx], length(idx), replace = FALSE)
    }
  }
  out
}

restricted_permutation_test <- function(data, response, predictor, strata, n_perm = 9999L) {
  model_data <- data %>%
    transmute(
      response = suppressWarnings(as.numeric(.data[[response]])),
      predictor = as.factor(.data[[predictor]]),
      strata = as.factor(.data[[strata]])
    ) %>%
    filter(is.finite(response), !is.na(predictor), !is.na(strata)) %>%
    droplevels()

  obs <- observed_f_stat(model_data, "response", "predictor", "strata")
  if (obs$status != "ok") {
    return(tibble(
      response = response,
      predictor = predictor,
      strata = strata,
      n_genomes = obs$n,
      n_strata = dplyr::n_distinct(model_data$strata),
      n_predictor_levels = dplyr::n_distinct(model_data$predictor),
      observed_F = obs$f_stat,
      df_num = obs$df_num,
      df_den = obs$df_den,
      p_value_permutation = NA_real_,
      n_perm = n_perm,
      status = obs$status
    ))
  }

  permuted_f <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    perm_data <- model_data
    perm_data$predictor <- permute_within_strata(perm_data$predictor, perm_data$strata)
    permuted_f[[i]] <- observed_f_stat(perm_data, "response", "predictor", "strata")$f_stat
  }

  valid_perm <- is.finite(permuted_f)
  p_perm <- (sum(permuted_f[valid_perm] >= obs$f_stat) + 1) / (sum(valid_perm) + 1)

  tibble(
    response = response,
    predictor = predictor,
    strata = strata,
    n_genomes = obs$n,
    n_strata = dplyr::n_distinct(model_data$strata),
    n_predictor_levels = dplyr::n_distinct(model_data$predictor),
    observed_F = obs$f_stat,
    df_num = obs$df_num,
    df_den = obs$df_den,
    p_value_permutation = p_perm,
    n_perm = n_perm,
    status = "ok"
  )
}

args <- parse_args()

genome_table <- get_arg(args, "genome-table", required = TRUE)
responses <- split_csv(get_arg(args, "responses", required = TRUE))
predictor <- get_arg(args, "predictor", "lifestyle")
strata <- get_arg(args, "strata", "phylum")
n_perm <- as_integer_arg(args, "n-perm", 9999L)
seed <- as_integer_arg(args, "seed", 1L)
output <- get_arg(args, "output", required = TRUE)
fdr_method <- get_arg(args, "fdr-method", "BH")

if (!file.exists(genome_table)) stop("Missing --genome-table: ", genome_table)
if (length(responses) == 0L) stop("No --responses supplied")
if (n_perm < 1L) stop("--n-perm must be >= 1")

set.seed(seed)

genomes <- read_tsv(genome_table, show_col_types = FALSE, progress = FALSE)
required_cols <- c(responses, predictor, strata)
missing_cols <- setdiff(required_cols, names(genomes))
if (length(missing_cols) > 0L) {
  stop("Genome table lacks columns: ", paste(missing_cols, collapse = ", "))
}

results <- bind_rows(lapply(responses, function(resp) {
  restricted_permutation_test(
    data = genomes,
    response = resp,
    predictor = predictor,
    strata = strata,
    n_perm = n_perm
  )
})) %>%
  mutate(
    q_value_permutation = p.adjust(p_value_permutation, method = fdr_method)
  )

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write_tsv(results, output, na = "NA")
message("Saved restricted permutation tests: ", output)
