############################################################
# Parallel tAI computation for fungal genomes
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
# nontrimmed with pseudo genes
path_to_cds    <- "./data/CDS/nontrimmed"
path_to_trna   <- "./data/tRNAscanse/old/"
path_to_output <- "./results/nontrimmed/tAI_per_genome.tsv"

#trimmed - no difference
# path_to_cds    <- "./data/CDS/trimmed"
# path_to_trna   <- "./data/tRNAscanse/old/"
# path_to_output <- "./results/trimmed/tAI_per_genome.tsv"


# nontrimmed without pseudo genes
path_to_cds    <- "./data/CDS/nontrimmed"
path_to_trna   <- "./data/tRNAscanse/nopseudo/"
path_to_output <- "./results/nopseudo/tAI_per_genome.tsv"



# Non-standard genetic codes
genetic_code_map <- c(
  "GCA_000182965.3" = "12"
)

default_genetic_code <- "1"

mapping_table <- read_tsv("./data/shortened_with_tAI.tsv",
                          show_col_types = FALSE)


############################
# Helper functions
############################

get_codon_table_for_assembly <- function(assemblyid) {
  code <- genetic_code_map[assemblyid]
  if (is.na(code)) code <- default_genetic_code
  get_codon_table(gcid = unname(code))
}

extract_assembly_id <- function(filename) {
  base <- basename(filename)
  id_with_version <- sub("cds_from_genomic.*", "", base)
  sub("\\.$", "", id_with_version)
}


read_trnascan_to_trna_level <- function(path) {

  df <- read_table2(
    path,
    skip = 3,
    col_names = c("seq", "trna_no", "begin", "end",
                  "type", "anticodon",
                  "intron_begin", "intron_end",
                  "score", "note"),
    col_types = cols(.default = col_character())
  )

  df_clean <- df %>%
    filter(is.na(note) | note == "") %>%
    mutate(trna_id = paste(type, anticodon, sep = "-")) %>%
    count(trna_id, name = "gcn")

  setNames(df_clean$gcn, df_clean$trna_id)
}


############################
# Core tAI computation
############################

compute_genome_tai <- function(fasta_file, trna_file, assemblyid) {

  codon_table <- get_codon_table_for_assembly(assemblyid)

  seqs <- readDNAStringSet(fasta_file)

  seq_qc <- check_cds(
    seqs,
    codon_table = codon_table,
    check_istop = TRUE,
    rm_start = TRUE,
    rm_stop  = TRUE
  )

  cf <- count_codons(seq_qc)

  trna_level <- read_trnascan_to_trna_level(trna_file)

  tRNA_w <- est_trna_weight(trna_level, codon_table)

  tai <- get_tai(cf, tRNA_w)

  mean(tai, na.rm = TRUE)
}



############################
# Parallel worker
############################

process_one_assembly <- function(assemblyid) {

  tryCatch({

    # Find CDS file
    fasta_file <- list.files(
      path = path_to_cds,
      pattern = paste0("^", assemblyid),
      full.names = TRUE
    )

    if (length(fasta_file) == 0) return(NULL)

    # Remove version suffix (.1, .2 etc.)
    assembly_short <- sub("\\..*", "", assemblyid)

    # Find corresponding tRNAscan file

    trna_file <- list.files(
      path = path_to_trna,
      pattern = paste0("^tRNAscan-SE_output_", assembly_short),
      full.names = TRUE
    )


    if (length(trna_file) == 0) {
      message("Missing tRNAscan: ", assemblyid)
      return(NULL)
    }

    mean_tai <- compute_genome_tai(
      fasta_file = fasta_file[1],
      trna_file  = trna_file[1],
      assemblyid = assemblyid
    )

    tibble(
      assemblyid = assemblyid,
      mean_tAI   = mean_tai
    )

  }, error = function(e) {
    message("FAILED: ", assemblyid, " | ", e$message)
    NULL
  })
}


############################
# Main execution
############################

n_cores <- max(1L, as.integer(detectCores() / 3))

results_list <- mclapply(
  mapping_table$assemblyid,
  process_one_assembly,
  mc.cores = n_cores
)

results_list <- Filter(Negate(is.null), results_list)

final_table <- bind_rows(results_list)

write_tsv(final_table, path_to_output)
