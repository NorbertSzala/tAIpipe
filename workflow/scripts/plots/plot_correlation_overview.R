suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

source(snakemake@input[["plotting_utils"]])
source(snakemake@input[["table_validation_utils"]])
source(snakemake@input[["plot_data_utils"]])
source(snakemake@input[["label_utils"]])

cfg <- read_plot_config(snakemake@input[["plotting_config"]])

plot_cor_matrix <- function(data, title) {
  require_columns(data, c("variable_x", "variable_y", "correlation"), "correlation table")

  ggplot(data, aes(x = variable_x, y = variable_y, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Spearman r",
      title = title
    ) +
    project_theme(cfg) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

gene <- read_tsv_checked(
  snakemake@input[["gene_correlations"]],
  required_columns = c("variable_x", "variable_y", "correlation"),
  table_name = "gene_level_correlations.tsv"
)

genome <- read_tsv_checked(
  snakemake@input[["genome_correlations"]],
  required_columns = c("variable_x", "variable_y", "correlation"),
  table_name = "genome_level_correlations.tsv"
)

gene_plot <- plot_cor_matrix(gene, "Gene-level correlations")
genome_plot <- plot_cor_matrix(genome, "Genome-level correlations")

save_plot(gene_plot, snakemake@output[["gene_png"]], cfg, size = "square")
save_plot(gene_plot, snakemake@output[["gene_pdf"]], cfg, size = "square")
save_plot(genome_plot, snakemake@output[["genome_png"]], cfg, size = "square")
save_plot(genome_plot, snakemake@output[["genome_pdf"]], cfg, size = "square")