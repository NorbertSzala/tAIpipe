# workflow/scripts/generate_report.R
log <- file(snakemake@log[[1]], open="wt")
sink(log); sink(log, type="message")

library(rmarkdown)

template_path <- "workflow/scripts/report_template.Rmd"

# Define the list of dynamic parameters directly from the snakemake object
report_params <- list(
  per_genome_dir  = normalizePath(snakemake@input$per_genome_dir),
  metadata_master = normalizePath(snakemake@input$metadata_master),
  samples_sheet   = normalizePath(snakemake@input$samples_sheet)
)

# 1. Render to static GitHub Flavored Markdown (.md)
rmarkdown::render(
  input         = template_path,
  output_format = "md_document",
  output_file   = basename(snakemake@output$md_report),
  output_dir    = dirname(snakemake@output$md_report),
  params        = report_params,
  quiet         = TRUE
)

# 2. Render to interactive HTML document
rmarkdown::render(
  input         = template_path,
  output_format = "html_document",
  output_file   = basename(snakemake@output$html_report),
  output_dir    = dirname(snakemake@output$html_report),
  params        = report_params,
  quiet         = TRUE
)

message("Report generation completed successfully.")