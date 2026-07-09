


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

source(snakemake@input[["plotting_utils"]])
source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")


cfg <- read_plot_config(snakemake@input[["plotting_config"]])

heatmap <- read_tsv_checked(
  snakemake@input[["heatmap"]],
  required_columns = c("sample", "codon", "value"),
  table_name = "codon_profile_heatmap.tsv"
)

variability <- read_tsv_checked(
  snakemake@input[["variability"]],
  required_columns = c("codon", "cv"),
  table_name = "codon_profile_variability.tsv"
)

heatmap_plot <- ggplot(heatmap, aes(x = codon, y = sample, fill = value)) +
  geom_tile() +
  scale_fill_continuous_project(cfg, option = "viridis") +
  labs(x = "Codon", y = "Genome", fill = "Value", title = "Codon profile heatmap") +
  project_theme(cfg) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

variability_plot <- variability |>
  arrange(desc(cv)) |>
  slice_head(n = 30) |>
  ggplot(aes(x = reorder(codon, cv), y = cv)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Coefficient of variation", title = "Most variable codons") +
  project_theme(cfg)

save_plot(heatmap_plot, snakemake@output[["trna_heatmap_png"]], cfg, size = "tall")
save_plot(heatmap_plot, snakemake@output[["trna_heatmap_pdf"]], cfg, size = "tall")
save_plot(variability_plot, snakemake@output[["variability_png"]], cfg, size = "wide")
save_plot(variability_plot, snakemake@output[["variability_pdf"]], cfg, size = "wide")
