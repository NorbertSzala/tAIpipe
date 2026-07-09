#!/usr/bin/env Rscript


# Builds the canonical long-format codon_profiles.tsv table. It combines codon counts, codon frequencies, RSCU values, tRNA weights, optimal-codon flags, genetic-code information, and sample metadata into one row per sample and codon.


suppressPackageStartupMessages({
    library(argparse)
    library(dplyr)
    library(purrr)
    library(readr)
    library(tidyr)
    library(tibble)
})



parser <- ArgumentParser(description = "Build canonical codon_profiles.tsv")
parser$add_argument("--metadata-dataset", required = TRUE)
parser$add_argument("--per-genome-dir", required = TRUE)
parser$add_argument("--output", required = TRUE)
args <- parser$parse_args()

as_flag <- function(x) {
    tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

read_required_csv <- function(path) {
    if (!file.exists(path)) stop("Missing codon-metric file: ", path)
    read_csv(path, show_col_types = FALSE, progress = FALSE)
}

prepare_rscu <- function(table, prefix) {
    required <- c("codon", "rscu", "w_cai")
    missing <- setdiff(required, names(table))
    if (length(missing) > 0L) {
        stop("RSCU table lacks columns: ", paste(missing, collapse = ", "))
    }

    metadata_columns <- intersect(c("codon", "aa_code", "amino_acid", "subfam"), names(table))
    table %>%
        select(all_of(metadata_columns), rscu, w_cai) %>%
        distinct(codon, .keep_all = TRUE) %>%
        rename(
            !!paste0(prefix, "_RSCU") := rscu,
            !!paste0(prefix, "_CAI_weight") := w_cai
        )
}

build_sample_profile <- function(metadata_row) {
    sample_id <- as.character(metadata_row$sample[[1]])
    metric_dir <- file.path(args$per_genome_dir, sample_id, "codon_metrics")

    counts <- read_required_csv(file.path(metric_dir, paste0(sample_id, "_codon_counts.csv")))
    genome_rscu <- read_required_csv(file.path(metric_dir, paste0(sample_id, "_rscu.csv")))
    reference_rscu <- read_required_csv(file.path(metric_dir, paste0(sample_id, "_reference_rscu.csv")))
    trna_weights <- read_required_csv(file.path(metric_dir, paste0(sample_id, "_trna_weights.csv")))

    if (!"seq_id" %in% names(counts)) {
        stop("Codon-count table lacks seq_id: ", sample_id)
    }
    codon_columns <- grep("^[ACGT]{3}$", names(counts), value = TRUE)
    if (length(codon_columns) == 0L) {
        stop("No DNA codon columns found in codon-count table: ", sample_id)
    }

    count_profile <- counts %>%
        summarise(across(all_of(codon_columns), ~ sum(as.numeric(.x), na.rm = TRUE))) %>%
        pivot_longer(everything(), names_to = "codon", values_to = "codon_count") %>%
        mutate(
            codon_count = as.numeric(codon_count),
            codon_frequency = if (sum(codon_count) > 0) codon_count / sum(codon_count) else NA_real_
        )

    genome_profile <- prepare_rscu(genome_rscu, "genome")
    reference_profile <- prepare_rscu(reference_rscu, "reference")

    required_trna <- c("codon", "anticodon", "trna_id", "ac_level", "W", "w")
    missing_trna <- setdiff(required_trna, names(trna_weights))
    if (length(missing_trna) > 0L) {
        stop("tRNA weight table lacks columns: ", paste(missing_trna, collapse = ", "))
    }
    trna_profile <- trna_weights %>%
        select(any_of(c(
            "codon", "anticodon", "trna_id", "ac_level", "W", "w",
            "aa_code", "amino_acid", "subfam"
        ))) %>%
        distinct(codon, .keep_all = TRUE) %>%
        rename(
            trna_anticodon = anticodon,
            trna_copy_number = ac_level,
            trna_absolute_weight = W,
            tRNA_weight = w
        )

    # Join annotation columns from the genome-wide RSCU table first. Duplicate
    # aa/subfamily columns from the other tables are excluded before joining.
    reference_values <- reference_profile %>%
        select(codon, reference_RSCU, reference_CAI_weight)
    trna_values <- trna_profile %>%
        select(codon, any_of(c(
            "trna_anticodon", "trna_id", "trna_copy_number",
            "trna_absolute_weight", "tRNA_weight"
        )))

    profile <- count_profile %>%
        left_join(genome_profile, by = "codon") %>%
        left_join(reference_values, by = "codon") %>%
        left_join(trna_values, by = "codon")

    # Stop codons may be present in raw count tables but do not have synonymous
    # codon/CAI/tRNA annotations. They are marked explicitly and retained.
    profile <- profile %>%
        mutate(is_stop_codon = codon %in% c("TAA", "TAG", "TGA")) %>%
        group_by(subfam) %>%
        mutate(
            reference_optimal_codon = if (
                !is.na(first(subfam)) && any(is.finite(reference_CAI_weight))
            ) {
                is.finite(reference_CAI_weight) &
                    reference_CAI_weight == max(reference_CAI_weight, na.rm = TRUE)
            } else {
                rep(FALSE, n())
            }
        ) %>%
        ungroup() %>%
        mutate(
            sample = sample_id,
            species = as.character(metadata_row$species[[1]]),
            accession = as.character(metadata_row$accession[[1]]),
            domain = as.character(metadata_row$domain[[1]]),
            kingdom = as.character(metadata_row$kingdom[[1]]),
            phylum = as.character(metadata_row$phylum[[1]]),
            lifestyle = as.character(metadata_row$lifestyle[[1]]),
            genetic_code = as.integer(metadata_row$genetic_code[[1]])
        ) %>%
        select(
            sample, species, accession, domain, kingdom, phylum, lifestyle,
            genetic_code, codon, is_stop_codon, aa_code, amino_acid, subfam,
            codon_count, codon_frequency, genome_RSCU, reference_RSCU,
            reference_CAI_weight, reference_optimal_codon,
            trna_anticodon, trna_id, trna_copy_number,
            trna_absolute_weight, tRNA_weight
        )

    if (anyDuplicated(profile$codon)) {
        stop("Duplicate codons in compiled profile for sample: ", sample_id)
    }
    profile
}

metadata <- read_tsv(args$metadata_dataset, show_col_types = FALSE, progress = FALSE)
required_metadata <- c(
    "sample", "species", "accession", "domain", "kingdom", "phylum",
    "lifestyle", "genetic_code", "include"
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

result <- map_dfr(seq_len(nrow(metadata)), ~ build_sample_profile(metadata[.x, , drop = FALSE]))
if (anyDuplicated(result[c("sample", "codon")])) {
    stop("Compiled codon profile contains duplicate sample + codon keys")
}

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
write_tsv(result, args$output, na = "NA")
message("Saved codon profile: ", args$output)
message("Rows: ", nrow(result))
