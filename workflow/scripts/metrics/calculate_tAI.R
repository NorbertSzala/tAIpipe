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
  library(cubar)
  library(Biostrings)
  library(readr)
  library(dplyr)
  library(optparse)
  library(tibble)
})


# -----------------
# --- Arguments ---
# -----------------
option_list <- list(
  make_option(
    c("-I", "--input"),
    type = "character",
    help = "input CDS fasta file"
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
    help = "Genetic code ID. Usually 1, C. albicans = 12 [default: %default]"
  ),
  make_option(
    c("-T", "--trna"),
    type = "character",
    help = "tRNA gene copy number file for tAI written as .tsv (converted tRNAscan-SE output)"
  ),
  make_option(
    c("-D", "--domain"),
    type = "character",
    default = "Eukarya",
    help = "Domain for tAI weights: Eukarya, Bacteria, or Archaea [default: %default]"
  ),
  make_option(
    c("-S", "--sample"),
    type = "character",
    help = "Sample name to produce varied outputs."
  )
)


parser <- OptionParser(
  usage = "%prog -I cds.fasta -O output_dir -G 1 -T trna_levels.tsv",
  option_list = option_list,
  description = "Calculate ENC, RSCU, CAI, tAI, AA usage, FOP, GC, and GC3s using cubar"
)

args <- parse_args(parser)


# ---------------------------
# --- Argument validation ---
# ---------------------------

if (is.null(args$input)) {
  stop("Missing input FASTA. Use -I / --input [path.fasta]")
}

if (!file.exists(args$input)) {
  stop("Input fasta does not exists: ", args$input)
}


if (is.null(args$trna)) {
  stop("Missing input tRNA gene copy number Use -T / --trna [path.tsv]")
}

if (!file.exists(args$trna)) {
  stop("tRNA file does not exist: ", args$trna)
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

if (!args$domain %in% c("Eukarya", "Bacteria", "Archaea")) {
  stop("--domain must be one of: Eukarya, Bacteria, Archaea")
}

if (is.null(args$sample)) {
  stop("Missing '--sample' or '-S' argument.")
}

if (is.null(args$genetic_code) || length(args$genetic_code) != 1 || is.na(args$genetic_code)) {
  stop("Invalid or missing genetic code. Use -G / --genetic-code, e.g. -G 1")
}


# -------------------------
# --- Helper functions  ---
# -------------------------

clean_seq_names <- function(x) {
  names(x) <- make.unique(names(x))
  x
}

load_cds <- function(path) {
  seq_raw <- readDNAStringSet(path)
  seq_raw <- clean_seq_names(seq_raw)

  return(seq_raw)
}

filter_valid_cds <- function(seq_raw, codon_table) {
  message("Input sequences: ", length(seq_raw))

  seq_checked <- check_cds(
    seq_raw,
    codon_table
  )

  message("Valid CDS after viltration: ", length(seq_checked))
  return(seq_checked)
}

write_table <- function(x, filename) {
  outpath <- file.path(args$outdir, filename)

  if (is.vector(x) && !is.list(x)) {
    df <- data.frame(
      seq_id = names(x),
      value = as.numeric(x),
      stringsAsFactors = FALSE
    )
  } else {
    df <- as.data.frame(x)
    df <- rownames_to_column(df, var = "seq_id")
  }

  readr::write_csv(df, outpath)
  message("Saved: ", outpath)
}


read_trna_levels <- function(path) {
  # Expected format: TSV:
  # AA-anticodon \t count
  # Ala-GCA \t 12

  trna_df <- readr::read_delim(
    path,
    delim = "\t", col_names = TRUE, show_col_types = FALSE
  )
  colnames(trna_df)[1:2] <- c("anticodon_id", "count")

  trna_level <- as.numeric(trna_df$count)
  names(trna_level) <- trna_df$anticodon_id

  return(trna_level)
}

# ----------------------
# --- Main execution ---
# ----------------------
sample <- args$sample

message("Loading codon table: ", as.character(args$genetic_code))
codon_table <- get_codon_table(gcid = as.character(args$genetic_code)) # 1 for almost _all organisms, 12 for C albicans. Reach from argument

message("Reading CDS fasta")
seq_raw <- load_cds(args$input)

message("Reading CDS for reference genes")
reference_seq_raw <- load_cds(args$reference_cds)

message("Checking CDS")
seq <- filter_valid_cds(
  seq_raw = seq_raw,
  codon_table = codon_table
)

message("Checking reference CDS")
reference_seq <- filter_valid_cds(
  seq_raw = reference_seq,
  codon_table = codon_table
)

message("Counting codons")
cf <- count_codons(seq)
write_table(cf, paste0(sample, "_codon_counts.csv"))

message("Counting reference codons")
reference_cf <- count_codons(reference_seq)
write_table(reference_cf, paste0(sample, "_reference_codon_counts.csv"))

# --- ENC ---

message("Calculating ENC...")
enc <- get_enc(
  cf = cf,
  codon_table = codon_table
)

write_table(enc, paste0(sample, "_enc.csv"))


# --- RSCU ---

message("Calculating RSCU...")
rscu <- est_rscu(
  cf = cf,
  codon_table = codon_table
)

write_table(rscu, paste0(sample, "_rscu.csv"))

message("Calculating reference RSCU...")
reference_rscu <- est_rscu(
  cf = reference_cf,
  codon_table = codon_table
)

write_table(reference_rscu, paste0(sample, "_reference_rscu.csv"))

# --- CAI ---

message("Calculating CAI...")
cai <- get_cai(
  cf = cf,
  rscu = rscu
)

write_table(cai, paste0(sample, "_cai.csv"))

message("Calculating reference CAI...")
reference_cai <- get_cai(
  cf = reference_cf,
  rscu = reference_rscu
)
write_table(reference_cai, paste0(sample, "_reference_cai.csv"))


# --- tAI ---

message("Reading tRNA levels...")
trna_level <- read_trna_levels(args$trna)

message("Estimating tRNA weights...")
trna_w <- est_trna_weight(
  trna_level = trna_level,
  codon_table = codon_table,
  domain = args$domain
)

# Controls
n_weights <- nrow(trna_w)
n_finite_weights <- sum(is.finite(trna_w$w))
n_positive_weights <- sum(is.finite(trna_w$w) & trna_w$w > 0)
if (n_positive_weights == 0L) {
  stop("No positive finite tRNA weights were estimated")
}

codon_totals <- colSums(cf)

valid_codons <- trna_w$codon[
  is.finite(trna_w$w) & trna_w$w > 0
]

used_codon_coverage <- (
  sum(codon_totals[names(codon_totals) %in% valid_codons]) /
    sum(codon_totals)
)

write_table(trna_w, paste0(sample, "_trna_weights.csv"))

message("Calculating tAI...")
tai <- get_tai(
  cf = cf,
  trna_w = trna_w
)

write_table(tai, paste0(sample, "_tai.csv"))


# --- Amino Acid Usage ---

message("Calculating amino acid usage...")
aau <- get_aau(
  cf = cf,
  codon_table = codon_table
)

write_table(aau, paste0(sample, "_amino_acid_usage.csv"))


# --- FOP ---

message("Calculating FOP...")
fop <- get_fop(
  cf = cf,
  codon_table = codon_table
)

write_table(fop, paste0(sample, "_fop.csv"))


# --- GC ---

message("Calculating GC...")
gc <- get_gc(cf = cf)

write_table(gc, paste0(sample, "_gc.csv"))


# --- GC3s ---

message("Calculating GC3s...")
gc3s <- get_gc3s(
  cf = cf,
  codon_table = codon_table
)

write_table(gc3s, paste0(sample, "_gc3s.csv"))


# ----------------------
# --- Summary table ---
# ----------------------

message("Creating summary...")

summary_df <- data.frame(
  seq_id = names(enc),
  ENC = as.numeric(enc),
  CAI = as.numeric(cai),
  FOP = as.numeric(fop),
  GC = as.numeric(gc),
  GC3s = as.numeric(gc3s),
  stringsAsFactors = FALSE
)

if (exists("tai")) {
  summary_df$tAI <- as.numeric(tai)
}

summary_path <- file.path(args$outdir, paste0(sample, "_summary.tsv"))
readr::write_tsv(summary_df, summary_path)

message("Saved: ", summary_path)
message("Done")
