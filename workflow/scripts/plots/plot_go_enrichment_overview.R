suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

source(snakemake@input[["plotting_utils"]])
source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")
source("workflow/scripts/lib/label_utils.R")



cfg <- read_plot_config(snakemake@input[["plotting_config"]])

terms <- read_tsv_checked(
  snakemake@input[["top_terms"]],
  table_name = "go_enrichment_top_terms.tsv"
)

p_col <- intersect(c("padj", "fdr", "q_value", "adjusted_p_value"), names(terms))[[1]]
effect_col <- intersect(c("odds_ratio", "estimate", "enrichment_ratio", "gene_ratio"), names(terms))[[1]]
term_col <- intersect(c("term_name", "go_name", "description", "go_term", "go_id"), names(terms))[[1]]

if (any(is.na(c(p_col, effect_col, term_col)))) {
  stop("GO plot requires term, effect, and adjusted p-value columns.", call. = FALSE)
}

dotplot <- terms %>%
  mutate(
    term_label = wrap_labels(.data[[term_col]], width = 45),
    minus_log10_padj = -log10(.data[[p_col]])
  ) %>%
  ggplot(aes(x = .data[[effect_col]], y = reorder(term_label, .data[[effect_col]]))) +
  geom_point(aes(size = minus_log10_padj)) +
  labs(x = effect_col, y = NULL, size = "-log10(FDR)", title = "Top GO enrichment terms") +
  project_theme(cfg)

forest <- terms %>%
  mutate(term_label = wrap_labels(.data[[term_col]], width = 45)) %>%
  ggplot(aes(x = .data[[effect_col]], y = reorder(term_label, .data[[effect_col]]))) +
  geom_vline(xintercept = 1, linewidth = 0.3) +
  geom_point() +
  labs(x = effect_col, y = NULL, title = "GO enrichment effect sizes") +
  project_theme(cfg)

save_plot(dotplot, snakemake@output[["dotplot_png"]], cfg, size = "tall")
save_plot(dotplot, snakemake@output[["dotplot_pdf"]], cfg, size = "tall")
save_plot(forest, snakemake@output[["forest_png"]], cfg, size = "tall")
save_plot(forest, snakemake@output[["forest_pdf"]], cfg, size = "tall")
