#!/usr/bin/env Rscript

# Reads CDS sequences, tRNA copy numbers, and the sample-specific genetic code, validates coding sequences, and calculates codon-usage metrics with the cubar package. It produces per-gene metrics and supporting codon-level outputs, and should use a validated ribosomal-gene reference set when estimating RSCU for CAI.


# Arguments:
# -I / --input        CDS FASTA
# -T / --trna         tRNA levels in  format:  AminoAcid-Anticodon e.g. Ala-GCA
# -O / --outdir       output directory
# -G / --genetic-code genetic code, np. 1 albo 12
# -D / --domain       Eukarya / Bacteria / Archaea


# Cubar package instruction:
# https://cran.r-project.org/web/packages/cubar/cubar.pdf


# -----------------
# --- Libraries ---
# -----------------

suppressPackageStartupMessages({
  library(Biostrings)
  library(cubar)
  library(dplyr)
  library(optparse)
  library(readr)
  library(tibble)
})


option_list <- list(
  make_option(c("-I", "--input"), type = "character", help = "Input CDS FASTA"),
  make_option(
    c("-R", "--reference-cds"),
    type = "character",
    dest = "reference_cds",
    help = "Ribosomal reference CDS FASTA used to estimate CAI weights"
  ),
  make_option(
    c("-O", "--outdir"),
    type = "character",
    default = "results/codon_usage_measurements",
    help = "Output directory [default: %default]"
  ),
  make_option(
    c("-G", "--genetic-code"),
    type = "integer",
    default = 1,
    dest = "genetic_code",
    help = "NCBI genetic-code ID [default: %default]"
  ),
  make_option(
    c("-T", "--trna"),
    type = "character",
    help = "Two-column tRNA copy-number TSV"
  ),
  make_option(
    c("-D", "--domain"),
    type = "character",
    default = "Eukarya",
    help = "Eukarya, Bacteria, or Archaea [default: %default]"
  ),
  make_option(c("-S", "--sample"), type = "character", help = "Sample ID"),
  make_option(
    "--min-reference-cds",
    type = "integer",
    default = 20L,
    dest = "min_reference_cds",
    help = "Warn when fewer valid reference CDS remain [default: %default]"
  ),
  make_option(
    "--qc-output",
    type = "character",
    dest = "qc_output",
    help = "Output TSV with per-sample metric QC"
  ),
  make_option(
    "--min-finite-tai-fraction",
    type = "double",
    default = 0.90,
    dest = "min_finite_tai_fraction"
  ),
  make_option(
    "--min-used-codon-coverage",
    type = "double",
    default = 0.95,
    dest = "min_used_codon_coverage"
  ),
  make_option(
    "--qc-mode",
    type = "character",
    default = "warn",
    dest = "qc_mode",
    help = "error, warn, or ignore"
  )
)

parser <- OptionParser(
  usage = paste(
    "%prog -I cds.fna -R ribosomal_reference.fna",
    "-T trna_counts.tsv -O output -S sample"
  ),
  option_list = option_list,
  description = "Calculate ENC, RSCU, CAI, tAI, amino-acid usage, FOP, GC and GC3s"
)
args <- parse_args(parser)

required_paths <- c(input = args$input, reference_cds = args$reference_cds, trna = args$trna)
for (label in names(required_paths)) {
  value <- required_paths[[label]]
  if (is.null(value) || is.na(value) || !nzchar(value)) {
    stop("Missing required argument: ", label)
  }
  if (!file.exists(value)) {
    stop("Input file does not exist (", label, "): ", value)
  }
}
if (is.null(args$sample) || is.na(args$sample) || !nzchar(args$sample)) {
  stop("Missing --sample argument")
}
if (!args$domain %in% c("Eukarya", "Bacteria", "Archaea")) {
  stop("--domain must be one of: Eukarya, Bacteria, Archaea")
}
if (is.null(args$genetic_code) || length(args$genetic_code) != 1L || is.na(args$genetic_code)) {
  stop("Invalid --genetic-code value")
}
if (is.null(args$min_reference_cds) || is.na(args$min_reference_cds) || args$min_reference_cds < 1L) {
  stop("--min-reference-cds must be a positive integer")
}
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

normalize_fasta_names <- function(sequences, path) {
  raw_names <- names(sequences)
  if (is.null(raw_names) || any(!nzchar(raw_names))) {
    stop("All FASTA records must have non-empty identifiers: ", path)
  }

  # cubar keeps full Biostrings record names. Duplicate full names would make
  # downstream joins ambiguous and are therefore rejected rather than renamed.
  duplicates <- unique(raw_names[duplicated(raw_names)])
  if (length(duplicates) > 0L) {
    stop(
      "Duplicated FASTA record names in ", path, ": ",
      paste(head(duplicates, 20L), collapse = ", ")
    )
  }
  sequences
}

load_cds <- function(path) {
  sequences <- readDNAStringSet(path)
  if (length(sequences) == 0L) {
    stop("No FASTA records found in: ", path)
  }
  normalize_fasta_names(sequences, path)
}

filter_valid_cds <- function(sequences, codon_table, label) {
  message(label, " input sequences: ", length(sequences))
  checked <- check_cds(sequences, codon_table)
  message(label, " valid CDS: ", length(checked))
  if (length(checked) == 0L) {
    stop("No valid CDS remained after cubar::check_cds for ", label)
  }
  checked
}

write_table <- function(x, filename) {
  outpath <- file.path(args$outdir, filename)

  if (is.vector(x) && !is.list(x)) {
    output <- tibble(
      seq_id = names(x),
      value = as.numeric(x)
    )
  } else {
    output <- as.data.frame(x, stringsAsFactors = FALSE)
    output <- rownames_to_column(output, var = "seq_id")
  }

  write_csv(output, outpath, na = "NA")
  message("Saved: ", outpath)
}

read_trna_levels <- function(path) {
  table <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  if (ncol(table) < 2L) {
    stop("tRNA table must contain at least two columns: identifier and count")
  }

  table <- table[, 1:2]
  names(table) <- c("anticodon_id", "count")
  table <- table %>%
    mutate(
      anticodon_id = trimws(as.character(anticodon_id)),
      count = suppressWarnings(as.numeric(count))
    )

  if (any(!nzchar(table$anticodon_id))) {
    stop("Empty anticodon identifiers found in: ", path)
  }
  if (any(!is.finite(table$count)) || any(table$count < 0)) {
    stop("tRNA counts must be finite, non-negative numbers: ", path)
  }
  if (anyDuplicated(table$anticodon_id)) {
    duplicates <- unique(table$anticodon_id[duplicated(table$anticodon_id)])
    stop("Duplicated tRNA identifiers: ", paste(head(duplicates, 20L), collapse = ", "))
  }

  levels <- table$count
  names(levels) <- table$anticodon_id
  levels
}

sample <- args$sample
message("Loading genetic code: ", args$genetic_code)
codon_table <- get_codon_table(gcid = as.character(args$genetic_code))

all_sequences_raw <- load_cds(args$input)
reference_sequences_raw <- load_cds(args$reference_cds)

all_sequences <- filter_valid_cds(
  all_sequences_raw, codon_table,
  label = "All CDS"
)

reference_sequences <- filter_valid_cds(
  reference_sequences_raw, codon_table,
  label = "Reference CDS"
)

if (length(reference_sequences) < args$min_reference_cds) {
  warning(
    "Only ", length(reference_sequences), " valid reference CDS remained; ",
    "configured warning threshold is ", args$min_reference_cds
  )
}
if (length(reference_sequences) < 5L) {
  stop("Fewer than five valid reference CDS are insufficient for stable CAI weights")
}

message("Counting codons")
cf <- count_codons(all_sequences)
reference_cf <- count_codons(reference_sequences)
write_table(cf, paste0(sample, "_codon_counts.csv"))

message("Calculating ENC")
enc <- get_enc(cf = cf, codon_table = codon_table)
write_table(enc, paste0(sample, "_enc.csv"))

message("Calculating genome-wide RSCU")
rscu <- est_rscu(cf = cf, codon_table = codon_table)
write_table(rscu, paste0(sample, "_rscu.csv"))

message("Estimating CAI weights from ribosomal reference CDS")
reference_rscu <- est_rscu(cf = reference_cf, codon_table = codon_table)
write_table(reference_rscu, paste0(sample, "_reference_rscu.csv"))

message("Calculating CAI with reference-derived weights")
cai <- get_cai(cf = cf, rscu = reference_rscu)
write_table(cai, paste0(sample, "_cai.csv"))

message("Reading tRNA levels")
trna_level <- read_trna_levels(args$trna)
message("Estimating tRNA weights")
trna_w <- est_trna_weight(
  trna_level = trna_level,
  codon_table = codon_table,
  domain = args$domain
)

if (!all(c("codon", "w") %in% names(trna_w))) {
  stop("cubar::est_trna_weight returned an unexpected table structure")
}
positive_weight <- is.finite(trna_w$w) & trna_w$w > 0
if (!any(positive_weight)) {
  stop("No positive finite tRNA weights were estimated")
}

codon_totals <- colSums(cf, na.rm = TRUE)
valid_codons <- unique(as.character(trna_w$codon[positive_weight]))
total_codons <- sum(codon_totals, na.rm = TRUE)
used_codon_coverage <- if (total_codons > 0) {
  sum(codon_totals[names(codon_totals) %in% valid_codons], na.rm = TRUE) / total_codons
} else {
  NA_real_
}
message("Used-codon coverage by positive tRNA weights: ", signif(used_codon_coverage, 5))
write_table(trna_w, paste0(sample, "_trna_weights.csv"))

message("Calculating tAI")
tai <- get_tai(cf = cf, trna_w = trna_w)
write_table(tai, paste0(sample, "_tai.csv"))

message("Calculating amino-acid usage")
aau <- get_aau(cf = cf, codon_table = codon_table)
write_table(aau, paste0(sample, "_amino_acid_usage.csv"))

message("Calculating FOP")
fop <- get_fop(cf = cf, codon_table = codon_table)
write_table(fop, paste0(sample, "_fop.csv"))

message("Calculating GC")
gc <- get_gc(cf = cf)
write_table(gc, paste0(sample, "_gc.csv"))

message("Calculating GC3s")
gc3s <- get_gc3s(cf = cf, codon_table = codon_table)
write_table(gc3s, paste0(sample, "_gc3s.csv"))

metric_vectors <- list(ENC = enc, CAI = cai, FOP = fop, GC = gc, GC3s = gc3s, tAI = tai)
expected_ids <- names(enc)
for (metric_name in names(metric_vectors)) {
  values <- metric_vectors[[metric_name]]
  if (!identical(names(values), expected_ids)) {
    stop("Identifier/order mismatch in metric: ", metric_name)
  }
}

summary_df <- tibble(
  seq_id = expected_ids,
  ENC = as.numeric(enc),
  CAI = as.numeric(cai),
  FOP = as.numeric(fop),
  GC = as.numeric(gc),
  GC3s = as.numeric(gc3s),
  tAI = as.numeric(tai)
)

finite_fraction <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L) return(NA_real_)
  mean(is.finite(x))
}

failed_checks <- character()

finite_tai_fraction <- finite_fraction(tai)
finite_cai_fraction <- finite_fraction(cai)
finite_enc_fraction <- finite_fraction(enc)

valid_cds_fraction <- length(all_sequences) / length(all_sequences_raw)
valid_reference_cds_fraction <- length(reference_sequences) / length(reference_sequences_raw)

if (!is.finite(finite_tai_fraction) ||
    finite_tai_fraction < args$min_finite_tai_fraction) {
  failed_checks <- c(
    failed_checks,
    paste0(
      "finite_tai_fraction<",
      args$min_finite_tai_fraction
    )
  )
}

if (!is.finite(used_codon_coverage) ||
    used_codon_coverage < args$min_used_codon_coverage) {
  failed_checks <- c(
    failed_checks,
    paste0(
      "used_codon_coverage<",
      args$min_used_codon_coverage
    )
  )
}

if (length(reference_sequences) < args$min_reference_cds) {
  failed_checks <- c(
    failed_checks,
    paste0("n_reference_cds<", args$min_reference_cds)
  )
}

qc_status <- if (length(failed_checks) == 0L) "PASS" else "WARN"

qc_df <- tibble(
  sample = sample,
  n_cds_input = length(all_sequences_raw),
  n_cds_valid = length(all_sequences),
  valid_cds_fraction = valid_cds_fraction,
  n_reference_cds_input = length(reference_sequences_raw),
  n_reference_cds_valid = length(reference_sequences),
  valid_reference_cds_fraction = valid_reference_cds_fraction,
  finite_tai_fraction = finite_tai_fraction,
  finite_cai_fraction = finite_cai_fraction,
  finite_enc_fraction = finite_enc_fraction,
  used_codon_coverage = used_codon_coverage,
  min_finite_tai_fraction = args$min_finite_tai_fraction,
  min_used_codon_coverage = args$min_used_codon_coverage,
  min_reference_cds = args$min_reference_cds,
  qc_status = qc_status,
  qc_reasons = if (length(failed_checks) == 0L) "none" else paste(failed_checks, collapse = ";")
)

if (!is.null(args$qc_output) && nzchar(args$qc_output)) {
  dir.create(dirname(args$qc_output), recursive = TRUE, showWarnings = FALSE)
  write_tsv(qc_df, args$qc_output, na = "NA")
  message("Saved metric QC: ", args$qc_output)
}

if (qc_status != "PASS" && args$qc_mode == "error") {
  stop("Metric QC failed for sample ", sample, ": ", qc_df$qc_reasons[[1]])
}


summary_path <- file.path(args$outdir, paste0(sample, "_summary.tsv"))
write_tsv(summary_df, summary_path, na = "NA")
message("Saved: ", summary_path)
message("Done")
