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

metrics <- read_tsv_checked(
  snakemake@input[["metric_summary"]],
  required_columns = c("group_variable", "group_value", "metric", "median"),
  table_name = "genome_metric_summary.tsv"
)

make_group_plot <- function(data, group_name, title) {
  data |>
    filter(group_variable == group_name) |>
    ggplot(aes(x = group_value, y = median)) +
    geom_col() +
    facet_wrap(vars(metric), scales = "free_y") +
    labs(x = NULL, y = "Median", title = title) +
    project_theme(cfg) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

phylum_plot <- make_group_plot(metrics, "phylum", "Genome metrics by phylum")
lifestyle_plot <- make_group_plot(metrics, "lifestyle", "Genome metrics by lifestyle")

save_plot(phylum_plot, snakemake@output[["phylum_png"]], cfg, size = "wide")
save_plot(phylum_plot, snakemake@output[["phylum_pdf"]], cfg, size = "wide")
save_plot(lifestyle_plot, snakemake@output[["lifestyle_png"]], cfg, size = "wide")
save_plot(lifestyle_plot, snakemake@output[["lifestyle_pdf"]], cfg, size = "wide")
