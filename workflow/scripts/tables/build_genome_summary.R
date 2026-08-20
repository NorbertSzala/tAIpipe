# Builds the canonical genome_summary.tsv table containing one row per genome. It aggregates gene-level codon metrics, genome FASTA statistics, tRNA quality-control results, annotation coverage, and biological metadata into a compact dataset for between-genome analyses.


suppressPackageStartupMessages({
    library(argparse)
    library(Biostrings)
    library(dplyr)
    library(purrr)
    library(readr)
    library(tibble)
})


parser <- ArgumentParser(description = "Build canonical genome_summary.tsv")
parser$add_argument("--gene-features", required = TRUE)
parser$add_argument("--metadata-dataset", required = TRUE)
parser$add_argument("--genome-dir", required = TRUE)
parser$add_argument("--per-genome-dir", required = TRUE)
parser$add_argument("--output", required = TRUE)
args <- parser$parse_args()

as_flag <- function(x) {
    tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

resolve_single_file <- function(directory, pattern) {
    matches <- Sys.glob(file.path(directory, pattern))
    if (length(matches) == 0L) {
        stop("No genome file found for pattern: ", file.path(directory, pattern))
    }
    if (length(matches) > 1L) {
        stop(
            "Multiple genome files found for pattern: ", file.path(directory, pattern),
            "\n", paste(matches, collapse = "\n")
        )
    }
    matches[[1]]
}

safe_mean <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
}

safe_median <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else median(x)
}

safe_sum <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else sum(x)
}

safe_fraction <- function(x) {
    if (length(x) == 0L) {
        return(NA_real_)
    }
    values <- as.logical(x)
    values <- values[!is.na(values)]
    if (length(values) == 0L) NA_real_ else mean(values)
}

read_one_row_qc <- function(path, prefix) {
    if (!file.exists(path)) {
        return(tibble())
    }
    table <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
    if (nrow(table) != 1L) {
        stop("Expected exactly one QC row in: ", path)
    }
    names(table) <- ifelse(
        names(table) == "sample",
        "sample",
        paste0(prefix, names(table))
    )
    table
}

summarize_genome_fasta <- function(sample_id, genome_path) {
    genome <- readDNAStringSet(genome_path)
    if (length(genome) == 0L) {
        stop("No sequences found in genome FASTA: ", genome_path)
    }

    counts <- colSums(
        letterFrequency(genome, letters = c("A", "C", "G", "T", "N"), as.prob = FALSE)
    )
    canonical <- sum(counts[c("A", "C", "G", "T")])
    total <- sum(width(genome))
    contig_lengths <- sort(as.numeric(width(genome)), decreasing = TRUE)
    l50 <- which(cumsum(contig_lengths) >= total / 2)[[1]]

    tibble(
        sample = sample_id,
        genome_file = normalizePath(genome_path, mustWork = TRUE),
        n_contigs = length(genome),
        genome_size_bp = as.numeric(total),
        n50_bp = contig_lengths[[l50]],
        l50 = l50,
        genome_gc = if (canonical > 0) {
            as.numeric((counts[["G"]] + counts[["C"]]) / canonical)
        } else {
            NA_real_
        },
        genome_n_fraction = if (total > 0) as.numeric(counts[["N"]] / total) else NA_real_
    )
}

metadata <- read_tsv(args$metadata_dataset, show_col_types = FALSE, progress = FALSE)
required_metadata <- c(
    "sample", "species", "accession", "domain", "kingdom", "phylum",
    "lifestyle", "genetic_code", "genome", "include"
)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0L) {
    stop("Metadata lacks columns: ", paste(missing_metadata, collapse = ", "))
}
metadata <- metadata %>%
    mutate(include = as_flag(include)) %>%
    filter(include)
if (nrow(metadata) == 0L) stop("No included samples in metadata table")
if (anyDuplicated(metadata$sample)) stop("Sample IDs must be unique in metadata")

genes <- read_tsv(args$gene_features, show_col_types = FALSE, progress = FALSE)
required_gene_columns <- c("sample", "seq_id", "tAI", "ENC", "CAI", "GC3s")
missing_gene <- setdiff(required_gene_columns, names(genes))
if (length(missing_gene) > 0L) {
    stop("Gene feature table lacks columns: ", paste(missing_gene, collapse = ", "))
}

# Add optional columns as NA so the aggregation remains stable when external
# annotations are not supplied.
for (column in c(
    "FOP", "GC", "delta_ENC", "metrics_available", "cds_qc_pass",
    "signal_peptide_present", "tm_present", "lcr_present", "pfam_present",
    "go_terms", "trna_qc_pass", "protein_length_aa"
)) {
    if (!column %in% names(genes)) genes[[column]] <- NA
}

gene_summary <- genes %>%
    group_by(sample) %>%
    summarise(
        n_genes = n(),
        n_proteins = sum(is.finite(suppressWarnings(as.numeric(protein_length_aa)))),
        proteome_length_aa = safe_sum(protein_length_aa),
        n_genes_with_metrics = sum(as.logical(metrics_available), na.rm = TRUE),
        fraction_genes_with_metrics = safe_fraction(metrics_available),
        fraction_cds_qc_pass = safe_fraction(cds_qc_pass),
        mean_tAI = safe_mean(tAI),
        median_tAI = safe_median(tAI),
        mean_ENC = safe_mean(ENC),
        median_ENC = safe_median(ENC),
        mean_CAI = safe_mean(CAI),
        median_CAI = safe_median(CAI),
        mean_FOP = safe_mean(FOP),
        median_FOP = safe_median(FOP),
        mean_GC = safe_mean(GC),
        median_GC = safe_median(GC),
        mean_GC3s = safe_mean(GC3s),
        median_GC3s = safe_median(GC3s),
        mean_delta_ENC = safe_mean(delta_ENC),
        median_delta_ENC = safe_median(delta_ENC),
        fraction_signal_peptide = safe_fraction(signal_peptide_present),
        fraction_tm = safe_fraction(tm_present),
        fraction_lcr = safe_fraction(lcr_present),
        fraction_pfam = safe_fraction(pfam_present),
        fraction_with_go = mean(!is.na(go_terms) & nzchar(trimws(as.character(go_terms)))),
        trna_qc_pass = if (all(is.na(trna_qc_pass))) NA else all(as.logical(trna_qc_pass), na.rm = TRUE),
        .groups = "drop"
    )

missing_samples <- setdiff(metadata$sample, gene_summary$sample)
if (length(missing_samples) > 0L) {
    stop("Gene table has no rows for included samples: ", paste(missing_samples, collapse = ", "))
}
extra_samples <- setdiff(gene_summary$sample, metadata$sample)
if (length(extra_samples) > 0L) {
    stop("Gene table contains samples not included in metadata: ", paste(extra_samples, collapse = ", "))
}

genome_stats <- map2_dfr(
    metadata$sample,
    metadata$genome,
    ~ summarize_genome_fasta(.x, resolve_single_file(args$genome_dir, .y))
)

qc_tables <- map_dfr(metadata$sample, function(sample_id) {
    trna_path <- file.path(
        args$per_genome_dir, sample_id, "qc", paste0(sample_id, "_trna_profile_qc.tsv")
    )
    kofam_path <- file.path(
        args$per_genome_dir, sample_id, "kofamscan", "ribosome_qc.tsv"
    )

    metric_path <- file.path(
        args$per_genome_dir, sample_id, "qc", paste0(sample_id, "_metric_qc.tsv")
    )

    trna_qc <- read_one_row_qc(trna_path, "trna_")
    kofam_qc <- read_one_row_qc(kofam_path, "kofam_")
    metric_qc <- read_one_row_qc(metric_path, "metric_")

    base <- tibble(sample = sample_id)
    if (nrow(trna_qc) == 1L) base <- left_join(base, trna_qc, by = "sample")
    if (nrow(kofam_qc) == 1L) base <- left_join(base, kofam_qc, by = "sample")
    if (nrow(metric_qc) == 1L) base <- left_join(base, metric_qc, by = "sample")
    base
})

result <- metadata %>%
    select(-include, -genome) %>%
    left_join(genome_stats, by = "sample") %>%
    left_join(gene_summary, by = "sample") %>%
    left_join(qc_tables, by = "sample")

if (nrow(result) != nrow(metadata) || anyDuplicated(result$sample)) {
    stop("Genome summary must contain exactly one row per included sample")
}

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
write_tsv(result, args$output, na = "NA")
message("Saved genome summary: ", args$output)
message("Rows: ", nrow(result))
