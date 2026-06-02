# workflow/scripts/generate_report.R
log <- file(snakemake@log[[1]], open="wt")
sink(log); sink(log, type="message")

library(rmarkdown)

#TODO: HArdcoded path, lepiej dac przez params w config oraz w R przez normalizePath()
template_path = "workflow/scripts/report_template.Rmd"

#TODO: nie tworzysz folderow w outpucie, moze sie wywalic (mozliwe ze wczesniej juz to zawarłes, nei patrzylem)
# Define the clean, simplified parameters for the Rmd document environment
# Maps the precise Snakemake input aliases straight to the report workspace
report_params <- list(
  per_genome_dir = normalizePath(snakemake@params$per_genome_dir),
  dataset        = normalizePath(snakemake@input$dataset),
  samples        = normalizePath(snakemake@input$samples),
  output_dir     = normalizePath(dirname(snakemake@output$html_report))
)

# 1. Render to a static GitHub Flavored Markdown (.md) document
rmarkdown::render(
  input         = template_path,
  output_format = "md_document",
  output_file   = basename(snakemake@output$md_report),
  output_dir    = dirname(snakemake@output$md_report),
  params        = report_params,
  quiet         = TRUE
)

# 2. Render to an interactive HTML document
rmarkdown::render(
  input         = template_path,
  output_format = "html_document",
  output_file   = basename(snakemake@output$html_report),
  output_dir    = dirname(snakemake@output$html_report),
  params        = report_params,
  quiet         = TRUE
)

message("Report generation completed successfully.")