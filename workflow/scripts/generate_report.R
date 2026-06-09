# workflow/scripts/generate_report.R
log <- file(snakemake@log[[1]], open="wt")
sink(log); sink(log, type="message")

on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

library(rmarkdown)

template_path <- normalizePath(
  snakemake@params$template_path,
  mustWork = TRUE
)


# Create output directories explicitly.
html_output_dir <- dirname(snakemake@output$html_report)
md_output_dir <- dirname(snakemake@output$md_report)

dir.create(html_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(md_output_dir, recursive = TRUE, showWarnings = FALSE)

# Optional GO dictionary configuration.
go_dictionary_enabled <- isTRUE(snakemake@params$go_dictionary_enabled)
go_dictionary_path <- snakemake@params$go_dictionary_path

if (is.null(go_dictionary_path) || !go_dictionary_enabled) {
  go_dictionary_path <- ""
}

# Define the clean, simplified parameters for the Rmd document environment
# Maps the precise Snakemake input aliases straight to the report workspace
report_params <- list(
  per_genome_dir = normalizePath(snakemake@params$per_genome_dir, mustWork = TRUE),
  dataset        = normalizePath(snakemake@input$dataset, mustWork = TRUE),
  samples        = normalizePath(snakemake@input$samples, mustWork = TRUE),
  output_dir     = normalizePath(html_output_dir, mustWork = TRUE),

  go_dictionary_enabled = go_dictionary_enabled,
  go_dictionary_path = normalizePath(go_dictionary_path, mustWork = FALSE)
)

message("Rendering report from template: ", template_path)
message("GO dictionary enabled: ", report_params$go_dictionary_enabled)
message("GO dictionary path: ", report_params$go_dictionary_path)

# 1. Render to a static GitHub Flavored Markdown (.md) document
rmarkdown::render(
  input         = template_path,
  output_format = "md_document",
  output_file   = basename(snakemake@output$md_report),
  output_dir    = md_output_dir,
  params        = report_params,
  quiet         = TRUE,
  envir         = new.env(parent = globalenv())
)

# 2. Render to an interactive HTML document
rmarkdown::render(
  input         = template_path,
  output_format = "html_document",
  output_file   = basename(snakemake@output$html_report),
  output_dir    = html_output_dir,
  params        = report_params,
  quiet         = TRUE,
  envir         = new.env(parent = globalenv())
)

message("Report generation completed successfully.")