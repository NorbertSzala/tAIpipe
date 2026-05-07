############################################################
# README
#
# This script computes genome-wide codon usage statistics
# for fungal genomes based on CDS FASTA files.
#
# For each assembly, it calculates:
#   - CAI (mean, median, SD, Q1, Q3)
#   - RSCU (mean, median, SD, Q1, Q3, max)
#   - ENC (mean, median, SD, Q1, Q3)
#   - GC content (mean, median, SD)
#   - GC3 content (mean, median, SD)
#
# Calculations are performed per gene and then aggregated
# to genome level. Processing is parallelized over genomes.
#
# Requirements:
#   - cubar
#   - Biostrings
#   - dplyr, readr
#
############################################################


############################
# Libraries
############################
library(cubar)
library(Biostrings)
library(readr)
library(dplyr)
library(parallel)


############################
# Configuration
############################
#nontrimmed
path_to_cds           <- "./data/CDS/nontrimmed"
path_to_mapping_table <- "./data/shortened_with_tAI.tsv"
path_to_save_results  <- "./results/nontrimmed/CAI.tsv"

#trimmed
# path_to_cds           <- "./data/CDS/trimmed"
# path_to_mapping_table <- "./data/shortened_with_tAI.tsv"
# path_to_save_results  <- "./results/trimmed/CAI.tsv"



# Non-standard genetic codes (NCBI IDs)
genetic_code_map <- c(
  "GCA_000182965.3" = "12"  # Candida albicans
)

default_genetic_code <- "1"


############################
# Input data
############################
mapping_table <- read_tsv(path_to_mapping_table, show_col_types = FALSE)


############################
# Helper functions
############################

# Return appropriate codon table for a given assembly
get_codon_table_for_assembly <- function(assemblyid) {
  code <- genetic_code_map[assemblyid]
  if (is.na(code)) code <- default_genetic_code
  get_codon_table(gcid = unname(code))
}


############################
# Core computation
############################

# Compute genome-wide codon usage statistics for a single FASTA
count_genome_codon_stats <- function(fasta_file, assemblyid) {

  codon_table <- get_codon_table_for_assembly(assemblyid)
  cds <- readDNAStringSet(fasta_file)

  # Count codon frequencies
  cf <- count_codons(cds)

  # RSCU (genome-wide reference)
  rscu <- est_rscu(
    cf = cf,
    codon_table = codon_table,
    level = "amino_acid"
  )

  # CAI per gene
  cai <- get_cai(
    cf = cf,
    rscu = rscu,
    level = "amino_acid"
  )

  # ENC per gene
  enc <- get_enc(
    cf = cf,
    codon_table = codon_table,
    level = "subfam"
  )

  # GC and GC3 per gene
  gc_vals  <- get_gc(cf)
  gc3_vals <- get_gc3s(cf, codon_table, level = "amino_acid")

  # Quantiles
  q_cai  <- quantile(cai,  probs = c(0.25, 0.75), na.rm = TRUE)
  q_rscu <- quantile(rscu$rscu, probs = c(0.25, 0.75), na.rm = TRUE)
  q_enc  <- quantile(enc,  probs = c(0.25, 0.75), na.rm = TRUE)

  # Return named list of scalar values
  list(
    # CAI
    mean_CAI = mean(cai, na.rm = TRUE),
    median_CAI = median(cai, na.rm = TRUE),
    sd_CAI   = sd(cai, na.rm = TRUE),
    Q1_CAI   = q_cai[[1]],
    Q3_CAI   = q_cai[[2]],

    # RSCU
    mean_RSCU = mean(rscu$rscu, na.rm = TRUE),
    median_RSCU = median(rscu$rscu, na.rm = TRUE),
    sd_RSCU   = sd(rscu$rscu, na.rm = TRUE),
    Q1_RSCU   = q_rscu[[1]],
    Q3_RSCU   = q_rscu[[2]],
    max_RSCU  = max(rscu$rscu, na.rm = TRUE),

    # ENC
    mean_ENC = mean(enc, na.rm = TRUE),
    median_ENC = median(enc, na.rm = TRUE),
    sd_ENC   = sd(enc, na.rm = TRUE),
    Q1_ENC   = q_enc[[1]],
    Q3_ENC   = q_enc[[2]],

    # GC
    mean_GC = mean(gc_vals, na.rm = TRUE),
    median_GC = median(gc_vals, na.rm = TRUE),
    sd_GC   = sd(gc_vals, na.rm = TRUE),

    # GC3
    mean_GC3 = mean(gc3_vals, na.rm = TRUE),
    median_GC3 = median(gc3_vals, na.rm = TRUE),
    sd_GC3   = sd(gc3_vals, na.rm = TRUE)
  )
}


############################
# Parallel worker
############################

process_one_assembly <- function(assemblyid) {
  tryCatch({

    fasta_file <- list.files(
      path = path_to_cds,
      pattern = paste0("^", assemblyid),
      full.names = TRUE
    )

    if (length(fasta_file) == 0) return(NULL)

    stats <- count_genome_codon_stats(
      fasta_file = fasta_file[1],
      assemblyid = assemblyid
    )

    # Enforce numeric scalars (no list-columns)
    stats <- lapply(stats, function(x) {
      if (length(x) == 0 || all(is.na(x))) NA_real_ else as.numeric(x)[1]
    })

    tibble(assemblyid = assemblyid, !!!stats)

  }, error = function(e) {
    message("FAILED: ", assemblyid, " | ", e$message)
    NULL
  })
}


############################
# Main execution
############################

n_cores <- max(1L, as.integer(detectCores() / 3))

tai_table <- mapping_table %>%
  select(assemblyid, tAI)


results_list <- mclapply(
  mapping_table$assemblyid,
  process_one_assembly,
  mc.cores = n_cores
)

results_list <- Filter(Negate(is.null), results_list)

# Combine per-assembly results
cai_rscu_table <- bind_rows(results_list)

# Add tAI from input table
final_results <- cai_rscu_table %>%
  left_join(tai_table, by = "assemblyid")

# Save results
write_tsv(final_results, path_to_save_results)

# (optional) full table with metadata
mapping_table_with_stats <- mapping_table %>%
  left_join(final_results, by = "assemblyid")