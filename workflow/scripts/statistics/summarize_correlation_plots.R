suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")

as_character_vector <- function(x, default = character()) {
  if (is.null(x)) return(default)
  if (length(x) == 0L) return(default)
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  x <- as.character(x)
  if (length(x) == 1L && grepl(",", x, fixed = TRUE)) {
    x <- strsplit(x, ",", fixed = TRUE)[[1]]
  }
  x <- trimws(x)
  unique(x[nzchar(x)])
}

validate_method <- function(method) {
  method <- as.character(method)[[1]]
  allowed <- c("spearman", "pearson", "kendall")
  if (!method %in% allowed) {
    stop("Unsupported correlation method: ", method, ". Allowed: ", paste(allowed, collapse = ", "))
  }
  method
}

select_numeric_variables <- function(data, requested_vars, table_name) {
  present <- intersect(requested_vars, names(data))
  missing <- setdiff(requested_vars, names(data))
  if (length(missing) > 0L) {
    message(table_name, ": skipped missing variables: ", paste(missing, collapse = ", "))
  }

  numeric_present <- present[vapply(data[present], function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    finite <- x_num[is.finite(x_num)]
    length(finite) >= 3L && length(unique(finite)) >= 2L
  }, logical(1))]

  non_numeric_or_constant <- setdiff(present, numeric_present)
  if (length(non_numeric_or_constant) > 0L) {
    message(table_name, ": skipped non-numeric/constant variables: ", paste(non_numeric_or_constant, collapse = ", "))
  }

  numeric_present
}

cor_to_long <- function(data, requested_vars, level, method) {
  vars <- select_numeric_variables(data, requested_vars, paste0(level, " table"))
  if (length(vars) < 2L) {
    return(tibble(
      level = character(),
      variable_x = character(),
      variable_y = character(),
      correlation = numeric(),
      n_pairwise = integer(),
      method = character()
    ))
  }

  numeric_data <- data |>
    transmute(across(all_of(vars), ~ suppressWarnings(as.numeric(.x))))

  cor_mat <- stats::cor(numeric_data, use = "pairwise.complete.obs", method = method)

  n_mat <- outer(vars, vars, Vectorize(function(a, b) {
    sum(is.finite(numeric_data[[a]]) & is.finite(numeric_data[[b]]))
  }))
  dimnames(n_mat) <- list(vars, vars)

  as_tibble(as.data.frame(as.table(cor_mat))) |>
    rename(variable_x = Var1, variable_y = Var2, correlation = Freq) |>
    mutate(
      n_pairwise = as.integer(n_mat[cbind(as.character(variable_x), as.character(variable_y))]),
      level = level,
      method = method,
      .before = 1
    ) |>
    arrange(level, variable_x, variable_y)
}

genes <- read_tsv_checked(snakemake@input[["gene_features"]], table_name = "gene_features.tsv")
genomes <- read_tsv_checked(snakemake@input[["genome_summary"]], table_name = "genome_summary.tsv")

gene_vars <- as_character_vector(
  snakemake@params[["gene_variables"]],
  default = c("tAI", "CAI", "GC3s", "ENC", "protein_length_aa")
)
genome_vars <- as_character_vector(
  snakemake@params[["genome_variables"]],
  default = c("mean_tAI", "median_tAI", "mean_GC3s")
)
method <- validate_method(snakemake@params[["method"]])

gene_cor <- cor_to_long(genes, gene_vars, "gene", method)
genome_cor <- cor_to_long(genomes, genome_vars, "genome", method)

write_tsv_safe(gene_cor, snakemake@output[["gene_correlations"]])
write_tsv_safe(genome_cor, snakemake@output[["genome_correlations"]])
