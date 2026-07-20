# Creates clustered, upper-left triangular correlation matrices for gene- and
# genome-level summaries, with variables grouped by biological interpretation.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

source(snakemake@input[["plotting_utils"]])
source(snakemake@input[["table_validation_utils"]])
source(snakemake@input[["plot_data_utils"]])
source(snakemake@input[["label_utils"]])
source("workflow/scripts/lib/plot_style_helpers.R")

cfg <- read_plot_config(snakemake@input[["plotting_config"]])

build_symmetric_correlation_matrix <- function(data) {
  vars <- sort(unique(c(as.character(data$variable_x), as.character(data$variable_y))))
  mat <- matrix(NA_real_, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))

  for (i in seq_len(nrow(data))) {
    x <- as.character(data$variable_x[[i]])
    y <- as.character(data$variable_y[[i]])
    value <- suppressWarnings(as.numeric(data$correlation[[i]]))
    if (!x %in% vars || !y %in% vars || !is.finite(value)) next
    mat[x, y] <- value
    if (!is.finite(mat[y, x])) mat[y, x] <- value
  }

  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      if (!is.finite(mat[i, j]) && is.finite(mat[j, i])) mat[i, j] <- mat[j, i]
    }
  }
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1
  mat
}

correlation_family <- function(variable) {
  dplyr::case_when(
    variable %in% c(
      "tAI", "tAI_z", "mean_tAI", "median_tAI", "CAI", "mean_CAI",
      "median_CAI", "ENC", "mean_ENC", "median_ENC", "delta_ENC",
      "mean_delta_ENC", "median_delta_ENC"
    ) ~ "Codon adaptation",
    variable %in% c(
      "GC", "mean_GC", "median_GC", "GC3s", "mean_GC3s",
      "median_GC3s", "genome_gc", "genome_RSCU", "RSCU", "rscu"
    ) ~ "Sequence composition",
    stringr::str_detect(variable, "protein_length|lcr_|tm_") ~ "Protein structure",
    stringr::str_detect(variable, "^fraction_") ~ "Annotation coverage",
    TRUE ~ "Other"
  )
}

order_correlation_variables <- function(mat) {
  variables <- rownames(mat)
  family_levels <- c(
    "Codon adaptation", "Sequence composition", "Protein structure",
    "Annotation coverage", "Other"
  )
  families <- correlation_family(variables)
  ordered <- unlist(lapply(family_levels, function(family) {
    members <- variables[families == family]
    if (length(members) <= 1L) return(members)
    cluster_correlation_variables(mat[members, members, drop = FALSE])
  }), use.names = FALSE)
  ordered[ordered %in% variables]
}

plot_cor_matrix <- function(data, title, subtitle) {
  require_columns(data, c("variable_x", "variable_y", "correlation"), "correlation table")
  data <- data %>%
    transmute(
      variable_x = as.character(variable_x),
      variable_y = as.character(variable_y),
      correlation = suppressWarnings(as.numeric(correlation))
    ) %>%
    filter(!is.na(variable_x), !is.na(variable_y), is.finite(correlation))

  mat <- build_symmetric_correlation_matrix(data)
  order_vars <- order_correlation_variables(mat)
  ordered <- mat[order_vars, order_vars, drop = FALSE]
  display_rows <- rev(order_vars)

  d <- tidyr::expand_grid(
    row_variable = display_rows,
    column_variable = order_vars
  ) %>%
    mutate(
      matrix_row_index = match(row_variable, order_vars),
      display_row_index = match(row_variable, display_rows),
      column_index = match(column_variable, order_vars),
      correlation = ordered[cbind(matrix_row_index, column_index)]
    ) %>%
    # Reverse the row order and retain the cells above the anti-diagonal. The
    # resulting matrix occupies the upper-left half of the square.
    filter(column_index <= length(order_vars) - display_row_index + 1L) %>%
    mutate(
      column_variable = factor(column_variable, levels = order_vars),
      row_variable = factor(row_variable, levels = rev(display_rows))
    )

  ggplot(d, aes(x = column_variable, y = row_variable, fill = correlation)) +
    geom_tile(colour = "white", linewidth = 0.55) +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3.7, fontface = "bold") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1), na.value = "grey92",
      guide = guide_colourbar(
        title.position = "top",
        barwidth = grid::unit(22, "cm"),
        barheight = grid::unit(0.62, "cm")
      )
    ) +
    scale_x_discrete(labels = correlation_plain_labels, drop = FALSE) +
    scale_y_discrete(labels = correlation_plain_labels, drop = FALSE) +
    coord_fixed() +
    labs(
      x = NULL,
      y = NULL,
      fill = math_labels("spearman_rho"),
      title = title,
      subtitle = subtitle
    ) +
    project_theme(cfg) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12.5),
      axis.text.y = element_text(size = 12.5),
      plot.title = element_text(size = 17, face = "bold"),
      plot.subtitle = element_text(size = 11.3),
      legend.position = "bottom",
      legend.justification = "center",
      panel.grid = element_blank()
    )
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

gene_plot <- plot_cor_matrix(
  gene,
  "Gene-level correlations across individual genes",
  "Each coefficient uses individual gene rows pooled across genomes; these variables are not genome means."
)
genome_plot <- plot_cor_matrix(
  genome,
  "Genome-level correlations of per-genome summaries",
  "Labels beginning with 'Mean' are per-genome means; each observation is one genome."
)

save_plot(gene_plot, snakemake@output[["gene_png"]], cfg, size = "correlation_large")
save_plot(gene_plot, snakemake@output[["gene_pdf"]], cfg, size = "correlation_large")
save_plot(genome_plot, snakemake@output[["genome_png"]], cfg, size = "correlation_large")
save_plot(genome_plot, snakemake@output[["genome_pdf"]], cfg, size = "correlation_large")
