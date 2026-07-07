# Renders the optional final analysis report from prepared tables, statistical outputs, and previously generated figures. The script should not recompute the main analyses; its role is to assemble results into a reproducible HTML or Markdown document.

dir.create(dirname(snakemake@log[[1]]), recursive = TRUE, showWarnings = FALSE)
log_connection <- file(snakemake@log[[1]], open = "wt")
sink(log_connection)
sink(log_connection, type = "message")
on.exit(
  {
    sink(type = "message")
    sink()
    close(log_connection)
  },
  add = TRUE
)

suppressPackageStartupMessages(library(rmarkdown))

dir.create(
  snakemake@params$report_file_path,
  recursive = TRUE,
  showWarnings = FALSE
)

template_path <- normalizePath(snakemake@input$template, mustWork = TRUE)
report_params <- list(
  per_genome_dir = normalizePath(snakemake@params$per_genome_dir, mustWork = TRUE),
  dataset = normalizePath(snakemake@input$dataset, mustWork = TRUE),
  samples = normalizePath(snakemake@input$samples, mustWork = TRUE),
  output_dir = normalizePath(
    dirname(snakemake@output$html_report),
    mustWork = FALSE
  )
)

render(
  input = template_path,
  output_format = "md_document",
  output_file = basename(snakemake@output$md_report),
  output_dir = dirname(snakemake@output$md_report),
  params = report_params,
  quiet = TRUE,
  envir = new.env(parent = globalenv())
)

render(
  input = template_path,
  output_format = "html_document",
  output_file = basename(snakemake@output$html_report),
  output_dir = dirname(snakemake@output$html_report),
  params = report_params,
  quiet = TRUE,
  envir = new.env(parent = globalenv())
)

message("Report generation completed successfully.")
