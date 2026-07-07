#!/usr/bin/env Rscript

# Queries Ensembl BioMart for the requested gene attributes and converts the response into a standardized annotation table. The script normalizes identifiers and column names so that the downloaded annotations can be joined reliably with gene-level pipeline outputs.

suppressPackageStartupMessages({
    library(argparse)
    library(biomaRt)
    library(dplyr)
    library(readr)
    library(stringr)
    library(purrr)
    library(tibble)
})

parser <- ArgumentParser(description = "Fetch cached BioMart annotations for one sample")
parser$add_argument("--input-gene-table", required = TRUE)
parser$add_argument("--sample", required = TRUE)
parser$add_argument("--host", required = TRUE)
parser$add_argument("--mart", required = TRUE)
parser$add_argument("--dataset", required = TRUE)
parser$add_argument("--id-filter", required = TRUE)
parser$add_argument("--id-attribute", required = TRUE)
parser$add_argument("--uniprot-attribute", default = "uniprotkb_all")
parser$add_argument("--go-attribute", default = "go_id")
parser$add_argument("--chunk-size", type = "integer", default = 500L)
parser$add_argument("--retries", type = "integer", default = 3L)
parser$add_argument("--strict", action = "store_true", default = FALSE)
parser$add_argument("--output", required = TRUE)
args <- parser$parse_args()

collapse_unique <- function(x) {
    x <- unique(x[!is.na(x) & nzchar(x)])
    if (length(x) == 0L) NA_character_ else paste(sort(x), collapse = ";")
}

normalize_go_terms <- function(x) {
    hits <- unlist(
        str_extract_all(as.character(x), regex("GO:[0-9]{7}", ignore_case = TRUE)),
        use.names = FALSE
    )
    collapse_unique(toupper(hits))
}

write_empty <- function(path, sample_id) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write_tsv(
        tibble(
            sample = character(),
            protein_id = character(),
            uniprot_id = character(),
            go_terms = character(),
            biomart_host = character(),
            biomart_mart = character(),
            biomart_dataset = character(),
            retrieved_at_utc = character()
        ),
        path
    )
    message("No annotations written for sample: ", sample_id)
}

gene_table <- read_tsv(args$input_gene_table, show_col_types = FALSE)
required_columns <- c("sample", "protein_id")
missing_columns <- setdiff(required_columns, names(gene_table))
if (length(missing_columns) > 0L) {
    stop("Input gene table lacks columns: ", paste(missing_columns, collapse = ", "))
}

ids <- gene_table %>%
    filter(sample == args$sample, !is.na(protein_id), protein_id != "") %>%
    pull(protein_id) %>%
    unique()

if (length(ids) == 0L) {
    write_empty(args$output, args$sample)
    quit(save = "no", status = 0L)
}

query_once <- function(id_chunk) {
    mart <- useMart(
        biomart = args$mart,
        host = args$host,
        dataset = args$dataset
    )

    available_attributes <- listAttributes(mart)$name
    available_filters <- listFilters(mart)$name
    requested_attributes <- c(
        args$id_attribute,
        args$uniprot_attribute,
        args$go_attribute
    )

    missing_attributes <- setdiff(requested_attributes, available_attributes)
    if (length(missing_attributes) > 0L) {
        stop(
            "Dataset lacks requested BioMart attributes: ",
            paste(missing_attributes, collapse = ", ")
        )
    }
    if (!args$id_filter %in% available_filters) {
        stop("Dataset lacks requested BioMart filter: ", args$id_filter)
    }

    getBM(
        attributes = requested_attributes,
        filters = args$id_filter,
        values = id_chunk,
        mart = mart
    )
}

query_with_retry <- function(id_chunk) {
    last_error <- NULL
    for (attempt in seq_len(args$retries)) {
        result <- tryCatch(
            query_once(id_chunk),
            error = function(error) {
                last_error <<- error
                NULL
            }
        )
        if (!is.null(result)) {
            return(result)
        }
        Sys.sleep(min(2^attempt, 30))
    }
    stop(last_error)
}

chunks <- split(ids, ceiling(seq_along(ids) / args$chunk_size))
raw_results <- tryCatch(
    map_dfr(chunks, query_with_retry),
    error = function(error) {
        if (isTRUE(args$strict)) {
            stop(error)
        }
        warning("BioMart failed for ", args$sample, ": ", conditionMessage(error))
        NULL
    }
)

if (is.null(raw_results) || nrow(raw_results) == 0L) {
    write_empty(args$output, args$sample)
    quit(save = "no", status = 0L)
}

names(raw_results)[names(raw_results) == args$id_attribute] <- "protein_id"
names(raw_results)[names(raw_results) == args$uniprot_attribute] <- "uniprot_raw"
names(raw_results)[names(raw_results) == args$go_attribute] <- "go_raw"

annotations <- raw_results %>%
    filter(!is.na(protein_id), protein_id != "") %>%
    group_by(protein_id) %>%
    summarise(
        uniprot_id = collapse_unique(uniprot_raw),
        go_terms = normalize_go_terms(go_raw),
        .groups = "drop"
    ) %>%
    mutate(
        sample = args$sample,
        biomart_host = args$host,
        biomart_mart = args$mart,
        biomart_dataset = args$dataset,
        retrieved_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ) %>%
    select(
        sample,
        protein_id,
        uniprot_id,
        go_terms,
        biomart_host,
        biomart_mart,
        biomart_dataset,
        retrieved_at_utc
    )

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
write_tsv(annotations, args$output, na = "NA")
message("Saved BioMart annotations: ", args$output)
