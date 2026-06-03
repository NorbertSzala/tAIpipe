# workflow/scripts/generate_report.R
log <- file(snakemake@log[[1]], open="wt")
sink(log); sink(log, type="message")

library(rmarkdown)

#TODO: HArdcoded path, lepiej dac przez params w config oraz w R przez normalizePath()
#DONE: Created  report_config | template_path param which indicates path
template_path = snakemake@params$template_path


#TODO: nie tworzysz folderow w outpucie, moze sie wywalic (mozliwe ze wczesniej juz to zawarłes, nei patrzylem)
#DONE: W raporcie są tworzone foldery, i powinno to zawsze działąc, ale w rasie czego ...
if (!dir.exists(snakemake@params$report_file_path)) {
  dir.create(snakemake@params$report_file_path, recursive = TRUE)
}

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