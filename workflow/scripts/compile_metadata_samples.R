#!/usr/bin/env Rscript

# workflow/scripts/compile_metadata_samples.R
# Script aggregating CDS sequence metadata, fetching GO terms and UniProt online via biomaRt.

suppressPackageStartupMessages({
  library(argparse)
  library(Biostrings)
  library(biomaRt)
  library(stringr)
  library(dplyr)
  library(readr)
  library(purrr)
})

# 1. Parsing input arguments
parser <- ArgumentParser(description = "Compilation of sample and gene metadata for tAIpipe")
parser$add_argument("--metadata_dataset", required = TRUE, help = "Path to dataset_test.tsv file")
parser$add_argument("--per-genome-dir", required = TRUE, help = "Root directory with per_genome results")
parser$add_argument("--cds-dir", required = TRUE, help = "Directory with input CDS files (*.fna)")
parser$add_argument("--gcode-dir", required = TRUE, help = "Directory with genetic codes (kept for compatibility)")
parser$add_argument("--output", required = TRUE, help = "Output path for the consolidated samples.tsv file")

args <- parser$parse_args()

# Initialization of connection with BioMart (Ensembl Fungi)
message("[-] Initializing online connection with biomaRt (Ensembl Fungi)...")
mart <- useMart(biomart = "fungi_mart", host = "https://fungi.ensembl.org")

# Function for fetching GO and UniProt ID per specific organism dataset
fetch_biomart_annotations_by_dataset <- function(protein_ids, dataset_name) {
  tryCatch({
    valid_ids <- protein_ids[!is.na(protein_ids) & protein_ids != ""]
    if (length(valid_ids) == 0) return(NULL)
    
    message(paste("[-] Connecting to dataset:", dataset_name, "for", length(valid_ids), "unique proteins..."))
    ensembl_mart <- useDataset(dataset_name, mart = mart)
    
    # Using 'uniprotkb_all' for maximum cross-reference retrieval alongside 'go_id'
    go_data <- getBM(
      attributes = c('protein_id', 'uniprotkb_all', 'go_id'),
      filters = 'protein_id',
      values = valid_ids,
      mart = ensembl_mart
    )
    
    message(paste("[+] BioMart returned", nrow(go_data), "raw annotation rows for", dataset_name))
    
    if (nrow(go_data) == 0) return(NULL)
    
    # Aggregation and formatting of results
    annot_mapped <- go_data %>%
      filter(!is.na(protein_id) & protein_id != "") %>%
      group_by(protein_id) %>%
      summarise(
        uniprot_fetched = paste(unique(uniprotkb_all[!is.na(uniprotkb_all) & uniprotkb_all != ""]), collapse = ";"),
        go_fetched = paste(unique(go_id[!is.na(go_id) & go_id != ""]), collapse = ";"),
        .groups = "drop"
      ) %>%
      mutate(
        uniprot_fetched = ifelse(uniprot_fetched == "", NA_character_, uniprot_fetched),
        go_fetched = ifelse(go_fetched == "", NA_character_, go_fetched)
      )
    
    return(annot_mapped)
  }, error = function(e) {
    message(paste("[!] Warning: Failed to fetch data via biomaRt for dataset", dataset_name, ":", e$message))
    return(NULL)
  })
}

# 2. Loading the active dataset sheet
dataset_df <- read_tsv(args$metadata_dataset, show_col_types = FALSE) %>%
  filter(include == TRUE | include == "True")

if (nrow(dataset_df) == 0) {
  stop("[!] Error: No samples with the include=True flag found in the dataset file.")
}

all_samples_metadata <- list()

# 3. Processing loop for each sample (Extracting structure from FASTA)
for (i in 1:nrow(dataset_df)) {
  current_sample <- dataset_df$sample[i]
  cds_pattern    <- dataset_df$cds[i]
  
  message(paste("[-] Processing sample:", current_sample))
  
  cds_files <- Sys.glob(file.path(args$cds_dir, cds_pattern))
  if (length(cds_files) == 0) {
    warning(paste("[!] Warning: No CDS file found for pattern:", cds_pattern))
    next
  }
  cds_file <- cds_files[1]
  
  fasta_headers <- names(readDNAStringSet(cds_file))
  message(paste("    -> Found", length(fasta_headers), "sequences in FASTA file:", basename(cds_file)))
  
  parsed_cds <- map_df(fasta_headers, function(header) {
    
    gene_id      <- str_match(header, "\\[locus_tag=([^\\]]+)\\]")[,2]
    protein_id   <- str_match(header, "\\[protein_id=([^\\]]+)\\]")[,2]
    protein_name <- str_match(header, "\\[protein=([^\\]]+)\\]")[,2]
    location     <- str_match(header, "\\[location=([^\\]]+)\\]")[,2]
    db_xref      <- str_match(header, "\\[db_xref=([^\\]]+)\\]")[,2]
    
    # Extraction directly from the header (if it exists)
    uniprot_id   <- str_match(db_xref, "UniProtKB(?:/Swiss-Prot|/TrEMBL)?:([A-Za-z0-9]+)")[,2]
    if(is.na(uniprot_id)) {
      uniprot_id <- str_match(db_xref, "GOA:([A-Za-z0-9]+)")[,2]
    }
    
    go_matches  <- str_match_all(db_xref, "GO:([0-9]+)")[[1]]
    go_fallback <- if(nrow(go_matches) > 0) paste(unique(go_matches[,2]), collapse = ";") else NA_character_

    coords <- as.numeric(str_extract_all(location, "[0-9]+")[[1]])
    length_nt <- NA
    length_aa <- NA
    
    if (length(coords) >= 2) {
      length_nt <- max(coords) - min(coords) + 1
      length_aa <- floor((length_nt - 3) / 3)
    }
    
    seq_id_metrics <- str_split(header, " ")[[1]][1]
    
    tibble(
      seq_id_metrics = seq_id_metrics,
      gene_id      = ifelse(is.na(gene_id), protein_id, gene_id),
      protein_id   = protein_id,
      uniprot_id   = uniprot_id,
      protein_name = protein_name,
      length_nt    = length_nt,
      length_aa    = length_aa,
      go_fallback  = go_fallback
    )
  })

  summary_file <- file.path(args$per_genome_dir, current_sample, "codon_metrics", paste0(current_sample, "_summary.tsv"))
  
  if (file.exists(summary_file)) {
    metrics_df <- read_tsv(summary_file, show_col_types = FALSE)
    metrics_df <- metrics_df %>%
      mutate(seq_id = str_trim(str_extract(seq_id, "^[^\\[]+")))

    message(paste("    -> Found", nrow(metrics_df), "metric rows inside summary tsv."))

    sample_consolidated <- parsed_cds %>%
      inner_join(metrics_df, by = c("seq_id_metrics" = "seq_id")) %>%
      mutate(sample = current_sample)
    
    message(paste("    -> Successfully joined", nrow(sample_consolidated), "rows for sample:", current_sample))
    
    # Collect initial header statistics for debugging
    headers_with_go <- sum(!is.na(sample_consolidated$go_fallback))
    headers_with_uniprot <- sum(!is.na(sample_consolidated$uniprot_id))
    message(paste("    -> FASTA internal metadata tracking: GO terms found in", headers_with_go, 
                  "sequences; UniProt IDs found in", headers_with_uniprot, "sequences."))
    
    all_samples_metadata[[current_sample]] <- sample_consolidated
  } else {
    warning(paste("[!] Warning: No metrics summary file found:", summary_file))
  }
}

# Combining data from all samples before the batch split BioMart query
if (length(all_samples_metadata) == 0) {
  stop("[!] Error: Failed to consolidate data for any of the samples.")
}

full_df <- bind_rows(all_samples_metadata)
message(paste("[-] Total rows combined across all samples before external BioMart annotation:", nrow(full_df)))

# 4. ORGANISM-SPECIFIC BATCH BIOMART QUERY
# Map each sample name dynamically to its prospective Ensembl Fungi dataset format (e.g., Spombe -> spombe_eg_gene)
full_df <- full_df %>%
  mutate(target_dataset = paste0(tolower(sample), "_eg_gene"))

unique_datasets <- unique(full_df$target_dataset)
biomart_results_list <- list()

for (ds in unique_datasets) {
  ds_protein_ids <- full_df %>% 
    filter(target_dataset == ds) %>% 
    pull(protein_id) %>% 
    unique()
  
  ds_annots <- fetch_biomart_annotations_by_dataset(ds_protein_ids, ds)
  if (!is.null(ds_annots)) {
    biomart_results_list[[ds]] <- ds_annots %>% mutate(target_dataset = ds)
  }
}

# Combine all retrieved online records
if (length(biomart_results_list) > 0) {
  combined_biomart_df <- bind_rows(biomart_results_list)
  message(paste("[+] Total unique proteins successfully annotated online via BioMart across all datasets:", 
                nrow(combined_biomart_df)))
  
  message("[-] Integrating online fetched data with the compiled metadata table...")
  full_df <- full_df %>%
    left_join(combined_biomart_df, by = c("protein_id", "target_dataset")) %>%
    mutate(
      uniprot_id = case_when(
        !is.na(uniprot_id) ~ uniprot_id,
        !is.na(uniprot_fetched) ~ uniprot_fetched,
        TRUE ~ NA_character_
      ),
      go_terms = case_when(
        !is.na(go_fallback) ~ go_fallback,
        !is.na(go_fetched) ~ go_fetched,
        TRUE ~ NA_character_
      )
    )
} else {
  message("[!] Warning: Online integration with BioMart was skipped or returned zero results. Falling back exclusively to FASTA headers.")
  full_df <- full_df %>%
    mutate(go_terms = go_fallback)
}

# Final data completeness verification printout
final_go_count <- sum(!is.na(full_df$go_terms))
final_uniprot_count <- sum(!is.na(full_df$uniprot_id))
message(paste("[+] Final Output Stats: Total rows:", nrow(full_df), 
              "| Valid GO entries:", final_go_count, 
              "| Valid UniProt entries:", final_uniprot_count))

# Cleaning and final column formatting according to the expected template
final_metadata_df <- full_df %>%
  select(sample, gene_id, protein_id, uniprot_id, protein_name, length_nt, length_aa, go_terms, 
         ENC, CAI, FOP, tAI, GC, GC3s)

# Replacing explicit structural NA values with string literals for output consistency
final_metadata_df[is.na(final_metadata_df)] <- NA

ensure_dir <- dirname(args$output)
if (!dir.exists(ensure_dir)) dir.create(ensure_dir, recursive = TRUE)

readr::write_tsv(final_metadata_df, args$output, na = "NA")
message(paste("[+] Success! Consolidated metadata safely saved to:", args$output))