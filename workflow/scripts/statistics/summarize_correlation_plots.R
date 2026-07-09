suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")


gene_vars <- intersect(
  c("tAI", "CAI", "GC3s", "ENC", "delta_ENC", "protein_length_aa", "log_protein_length_aa"),
  names(readr::read_tsv(snakemake@input[["gene_features"]], n_max = 1, show_col_types = FALSE))
)

genome_vars <- intersect(
  c("mean_tAI", "median_tAI", "mean_GC3s", "median_delta_ENC", "n_genes", "n_trnas"),
  names(readr::read_tsv(snakemake@input[["genome_summary"]], n_max = 1, show_col_types = FALSE))
)

cor_to_long <- function(data, vars, level) {
  if (length(vars) < 2) {
    return(tibble::tibble())
  }

  mat <- stats::cor(data[vars], use = "pairwise.complete.obs", method = "spearman")

  as.data.frame(as.table(mat)) |>
    rename(variable_x = Var1, variable_y = Var2, correlation = Freq) |>
    mutate(level = level)
}

genes <- read_tsv_checked(snakemake@input[["gene_features"]], table_name = "gene_features.tsv")
genomes <- read_tsv_checked(snakemake@input[["genome_summary"]], table_name = "genome_summary.tsv")

gene_cor <- cor_to_long(genes, gene_vars, "gene")
genome_cor <- cor_to_long(genomes, genome_vars, "genome")

write_tsv_safe(gene_cor, snakemake@output[["gene_correlations"]])
write_tsv_safe(genome_cor, snakemake@output[["genome_correlations"]])
