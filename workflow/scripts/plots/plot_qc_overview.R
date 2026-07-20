suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source(snakemake@input[["plotting_utils"]])
source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_style_helpers.R")


cfg <- read_plot_config(snakemake@input[["plotting_config"]])

qc <- readr::read_tsv(snakemake@input[["qc_summary"]], show_col_types = FALSE)
require_non_empty(qc, "trna_qc_summary.tsv")

flag_cols <- grep("status|flag|result|qc", names(qc), ignore.case = TRUE, value = TRUE)

if (length(flag_cols) == 0) {
  stop("No QC flag/status columns found in QC summary.", call. = FALSE)
}

qc_long <- qc |>
  select(any_of(c("sample", flag_cols))) |>
  pivot_longer(cols = all_of(flag_cols), names_to = "qc_metric", values_to = "status")

qc_plot <- ggplot(qc_long, aes(x = qc_metric, fill = as.character(status))) +
  geom_bar(position = "stack") +
  coord_flip() +
  labs(x = NULL, y = "Count", fill = "Status", title = tex_label("$\\mathrm{tRNA}$ QC flags")) +
  project_theme(cfg)

save_plot(qc_plot, snakemake@output[["qc_png"]], cfg, size = "wide")
save_plot(qc_plot, snakemake@output[["qc_pdf"]], cfg, size = "wide")
