suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
})

source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")


qc_files <- unlist(snakemake@input[["qc_files"]], use.names = FALSE)

if (length(qc_files) == 0) {
  stop("No QC files were provided.", call. = FALSE)
}

qc <- purrr::map_dfr(
  qc_files,
  function(path) {
    data <- readr::read_tsv(path, show_col_types = FALSE)
    if (!"sample" %in% names(data)) {
      data$sample <- basename(dirname(dirname(path)))
    }
    data
  }
)

require_non_empty(qc, "merged QC table")

write_tsv_safe(qc, snakemake@output[["qc_summary"]])
