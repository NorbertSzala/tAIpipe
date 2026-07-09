#!/usr/bin/env Rscript

# Builds the canonical gene_features.tsv table containing one row per gene. It merges CDS metadata, per-gene codon metrics, sample metadata, tRNA quality-control status, and optional external annotations, and derives within-genome variables such as standardized tAI and tAI percentiles.


suppressPackageStartupMessages({
  library(argparse)
  library(Biostrings)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})



parser <- ArgumentParser(
  description = "Build the canonical gene-level feature table for tAIpipe"
)
parser$add_argument("--metadata-dataset", required = TRUE)
parser$add_argument("--per-genome-dir", required = TRUE)
parser$add_argument("--cds-dir", required = TRUE)
parser$add_argument(
  "--annotation-table",
  default = NULL,
  help = paste(
    "Optional TSV with a unique sample + protein_id or sample + seq_id key.",
    "All non-key annotation columns are preserved."
  )
)
parser$add_argument("--output", required = TRUE)
parser$add_argument(
  "--require-go-terms",
  action = "store_true",
  help = "Fail if no valid GO identifiers are present in the final gene table."
)

parser$add_argument(
  "--min-go-annotated-genes",
  type = "integer",
  default = 1L,
  help = "Minimum number of genes with valid GO terms when --require-go-terms is used."
)
args <- parser$parse_args()

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

resolve_single_file <- function(directory, pattern) {
  matches <- Sys.glob(file.path(directory, pattern))
  if (length(matches) == 0L) {
    stop("No file found for pattern: ", file.path(directory, pattern))
  }
  if (length(matches) > 1L) {
    stop(
      "Multiple files found for pattern: ", file.path(directory, pattern),
      "\n", paste(matches, collapse = "\n")
    )
  }
  matches[[1]]
}

normalize_seq_id <- function(x) {
  str_extract(str_trim(as.character(x)), "^[^[:space:]\\[]+")
}

extract_first_tag <- function(header, tag) {
  pattern <- paste0("\\[", tag, "=([^\\]]+)\\]")
  hit <- str_match(header, pattern)[, 2]
  ifelse(is.na(hit) | hit == "", NA_character_, hit)
}

collapse_unique <- function(x) {
  values <- unique(trimws(as.character(x)))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) NA_character_ else paste(sort(values), collapse = ";")
}

normalize_go_terms <- function(x) {
  if (length(x) == 0L) {
    return(NA_character_)
  }

  values <- as.character(x)
  values <- values[!is.na(values) & nzchar(trimws(values))]

  if (length(values) == 0L) {
    return(NA_character_)
  }

  hits <- unlist(
    str_extract_all(
      values,
      regex("GO:[0-9]{7}", ignore_case = TRUE)
    ),
    use.names = FALSE
  )

  hits <- toupper(hits)
  hits <- unique(hits[!is.na(hits) & nzchar(hits)])

  if (length(hits) == 0L) {
    return(NA_character_)
  }

  paste(sort(hits), collapse = ";")
}

extract_uniprot_ids <- function(header) {
  matches <- str_match_all(
    header,
    regex(
      "(?:UniProtKB(?:/Swiss-Prot|/TrEMBL)?|GOA):([A-Z0-9]+(?:-[0-9]+)?)",
      ignore_case = TRUE
    )
  )[[1]]

  if (nrow(matches) == 0L) {
    return(NA_character_)
  }
  collapse_unique(toupper(matches[, 2]))
}

merge_semicolon_values <- function(x, y) {
  values <- c(x, y)
  split_values <- unlist(
    str_split(values[!is.na(values)], ";", simplify = FALSE),
    use.names = FALSE
  )
  collapse_unique(split_values)
}

terminal_stop_status <- function(sequences, genetic_code_id) {
  widths <- width(sequences)
  last_codons <- rep(NA_character_, length(sequences))
  eligible <- widths >= 3L

  if (any(eligible)) {
    last_codons[eligible] <- as.character(
      subseq(
        sequences[eligible],
        start = widths[eligible] - 2L,
        end = widths[eligible]
      )
    )
  }

  genetic_code <- getGeneticCode(as.character(genetic_code_id))
  stop_codons <- names(genetic_code)[genetic_code == "*"]
  last_codons %in% stop_codons
}

safe_zscore <- function(x) {
  finite <- is.finite(x)
  result <- rep(NA_real_, length(x))

  if (sum(finite) < 2L) {
    return(result)
  }

  spread <- sd(x[finite])
  if (!is.finite(spread) || spread == 0) {
    return(result)
  }

  result[finite] <- (x[finite] - mean(x[finite])) / spread
  result
}

safe_percentile <- function(x) {
  finite <- is.finite(x)
  result <- rep(NA_real_, length(x))
  n <- sum(finite)

  if (n == 0L) {
    return(result)
  }
  if (n == 1L) {
    result[finite] <- 0.5
    return(result)
  }

  result[finite] <- (rank(x[finite], ties.method = "average") - 1) / (n - 1)
  result
}

parse_cds_fasta <- function(cds_file, genetic_code_id) {
  cds <- readDNAStringSet(cds_file)
  headers <- names(cds)

  if (is.null(headers) || any(!nzchar(headers))) {
    stop("All CDS FASTA records must have non-empty headers: ", cds_file)
  }

  seq_ids <- normalize_seq_id(headers)
  duplicated_ids <- unique(seq_ids[duplicated(seq_ids)])
  if (length(duplicated_ids) > 0L) {
    stop(
      "Duplicated normalized CDS identifiers in ", cds_file, ": ",
      paste(duplicated_ids, collapse = ", ")
    )
  }

  lengths_nt <- width(cds)
  has_terminal_stop <- terminal_stop_status(cds, genetic_code_id)
  lengths_aa <- pmax(
    0L,
    floor(lengths_nt / 3L) - ifelse(has_terminal_stop, 1L, 0L)
  )

  tibble(
    seq_id = seq_ids,
    gene_id_locus = map_chr(headers, extract_first_tag, tag = "locus_tag"),
    gene_id_gene = map_chr(headers, extract_first_tag, tag = "gene"),
    protein_id = map_chr(headers, extract_first_tag, tag = "protein_id"),
    protein_name = map_chr(headers, extract_first_tag, tag = "protein"),
    uniprot_id_header = map_chr(headers, extract_uniprot_ids),
    go_terms_header = map_chr(headers, normalize_go_terms),
    cds_length_nt = as.integer(lengths_nt),
    protein_length_aa = as.integer(lengths_aa),
    cds_length_multiple_of_three = lengths_nt %% 3L == 0L,
    has_terminal_stop = as.logical(has_terminal_stop)
  ) %>%
    mutate(
      gene_id = coalesce(gene_id_locus, gene_id_gene, protein_id, seq_id)
    ) %>%
    select(
      seq_id,
      gene_id,
      protein_id,
      protein_name,
      uniprot_id_header,
      go_terms_header,
      cds_length_nt,
      protein_length_aa,
      cds_length_multiple_of_three,
      has_terminal_stop
    )
}

read_metric_summary <- function(path) {
  metrics <- read_tsv(path, show_col_types = FALSE)
  if (!"seq_id" %in% names(metrics)) {
    stop("Metric summary lacks seq_id column: ", path)
  }

  metrics <- metrics %>%
    mutate(
      seq_id = normalize_seq_id(seq_id),
      metric_row_present = TRUE
    )

  duplicate_ids <- unique(metrics$seq_id[duplicated(metrics$seq_id)])
  if (length(duplicate_ids) > 0L) {
    stop(
      "Duplicated normalized sequence IDs in metric summary ", path, ": ",
      paste(duplicate_ids, collapse = ", ")
    )
  }
  metrics
}

read_trna_qc <- function(per_genome_dir, sample_id) {
  path <- file.path(
    per_genome_dir,
    sample_id,
    "qc",
    paste0(sample_id, "_trna_profile_qc.tsv")
  )

  if (!file.exists(path)) {
    stop("Missing tRNA QC table: ", path)
  }

  qc <- read_tsv(path, show_col_types = FALSE)
  if (nrow(qc) != 1L) {
    stop("Expected exactly one row in tRNA QC table: ", path)
  }

  required <- c("sample", "qc_status", "qc_reasons")
  missing <- setdiff(required, names(qc))
  if (length(missing) > 0L) {
    stop("tRNA QC table lacks columns: ", paste(missing, collapse = ", "))
  }

  qc %>%
    transmute(
      sample = as.character(sample),
      trna_qc_status = as.character(qc_status),
      trna_qc_reasons = as.character(qc_reasons),
      trna_qc_pass = qc_status == "PASS",
      n_total_trnas = as.integer(.data$n_total_trnas),
      n_elongator_trnas = as.integer(.data$n_elongator_trnas),
      n_unique_anticodons = as.integer(.data$n_unique_anticodons),
      n_trna_amino_acids = as.integer(.data$n_amino_acids)
    )
}

read_optional_annotations <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop("Annotation table does not exist: ", path)
  }

  annotations <- read_tsv(path, show_col_types = FALSE)
  if (!"sample" %in% names(annotations)) {
    stop("Annotation table must contain sample column: ", path)
  }

  annotation_payload_columns <- intersect(
    names(annotations),
    c(
      "go_terms",
      "uniprot_id",
      "pfam_terms",
      "signal_peptide_present",
      "tm_present",
      "lcr_present",
      "annotation_source"
    )
  )

  if (length(annotation_payload_columns) == 0L) {
    stop(
      "Annotation table contains join keys but no supported annotation columns. ",
      "Expected at least one of: go_terms, uniprot_id, pfam_terms, ",
      "signal_peptide_present, tm_present, lcr_present, annotation_source."
    )
  }

  if ("protein_id" %in% names(annotations)) {
    join_keys <- c("sample", "protein_id")
  } else if ("seq_id" %in% names(annotations)) {
    annotations <- annotations %>% mutate(seq_id = normalize_seq_id(seq_id))
    join_keys <- c("sample", "seq_id")
  } else {
    stop("Annotation table must contain protein_id or seq_id: ", path)
  }

  duplicate_keys <- annotations %>%
    count(across(all_of(join_keys)), name = "n") %>%
    filter(n > 1L)
  if (nrow(duplicate_keys) > 0L) {
    stop(
      "Annotation table must have unique ",
      paste(join_keys, collapse = " + "),
      " keys. Aggregate duplicate annotations before this step."
    )
  }

  if ("uniprot_id" %in% names(annotations)) {
    annotations <- annotations %>%
      rename(uniprot_id_external = uniprot_id)
  } else {
    annotations$uniprot_id_external <- NA_character_
  }

  if ("go_terms" %in% names(annotations)) {
    annotations <- annotations %>%
      mutate(go_terms_external = map_chr(go_terms, normalize_go_terms)) %>%
      select(-go_terms)

    n_go <- sum(!is.na(annotations$go_terms_external))
    if (n_go == 0L) {
      warning(
        "Annotation table has a go_terms column, but no valid GO identifiers ",
        "matching GO:[0-9]{7} were found."
      )
    }
  } else {
    annotations$go_terms_external <- NA_character_
  }

  list(data = annotations, keys = join_keys)
}



required_dataset_columns <- c(
  "sample", "species", "accession", "domain", "kingdom", "phylum",
  "lifestyle", "genetic_code", "cds", "include"
)

dataset <- read_tsv(args$metadata_dataset, show_col_types = FALSE)
missing_dataset_columns <- setdiff(required_dataset_columns, names(dataset))
if (length(missing_dataset_columns) > 0L) {
  stop(
    "Metadata dataset lacks columns: ",
    paste(missing_dataset_columns, collapse = ", ")
  )
}

dataset <- dataset %>%
  mutate(include = as_flag(include)) %>%
  filter(include)

if (nrow(dataset) == 0L) {
  stop("No samples with include=true in metadata dataset.")
}

external_annotations <- read_optional_annotations(args$annotation_table)
compiled <- vector("list", nrow(dataset))

for (i in seq_len(nrow(dataset))) {
  metadata_row <- dataset[i, , drop = FALSE]
  sample_id <- metadata_row$sample[[1]]
  genetic_code_id <- as.integer(metadata_row$genetic_code[[1]])
  cds_file <- resolve_single_file(args$cds_dir, metadata_row$cds[[1]])
  summary_file <- file.path(
    args$per_genome_dir,
    sample_id,
    "codon_metrics",
    paste0(sample_id, "_summary.tsv")
  )

  if (!file.exists(summary_file)) {
    stop("Missing codon metric summary: ", summary_file)
  }

  message("Processing sample: ", sample_id)

  local_metadata <- parse_cds_fasta(cds_file, genetic_code_id)
  metrics <- read_metric_summary(summary_file)
  trna_qc <- read_trna_qc(args$per_genome_dir, sample_id)

  sample_table <- local_metadata %>%
    left_join(metrics, by = "seq_id") %>%
    mutate(
      sample = sample_id,
      cds_qc_pass = !is.na(metric_row_present),
      metrics_available = if_any(
        any_of(c("ENC", "CAI", "FOP", "tAI", "GC", "GC3s")),
        ~ !is.na(.x)
      )
    ) %>%
    select(-metric_row_present) %>%
    left_join(trna_qc, by = "sample") %>%
    mutate(
      species = as.character(metadata_row$species[[1]]),
      accession = as.character(metadata_row$accession[[1]]),
      domain = as.character(metadata_row$domain[[1]]),
      kingdom = as.character(metadata_row$kingdom[[1]]),
      phylum = as.character(metadata_row$phylum[[1]]),
      lifestyle = as.character(metadata_row$lifestyle[[1]]),
      genetic_code = genetic_code_id
    )

  if (!is.null(external_annotations)) {
    sample_annotations <- external_annotations$data %>%
      filter(sample == sample_id)

    sample_table <- sample_table %>%
      left_join(sample_annotations, by = external_annotations$keys)
  } else {
    sample_table <- sample_table %>%
      mutate(
        uniprot_id_external = NA_character_,
        go_terms_external = NA_character_
      )
  }

  if (!"uniprot_id_external" %in% names(sample_table)) {
    sample_table$uniprot_id_external <- NA_character_
  }
  if (!"go_terms_external" %in% names(sample_table)) {
    sample_table$go_terms_external <- NA_character_
  }

  sample_table <- sample_table %>%
    rowwise() %>%
    mutate(
      uniprot_id = merge_semicolon_values(
        uniprot_id_header,
        uniprot_id_external
      ),
      go_terms = normalize_go_terms(
        c(go_terms_header, go_terms_external)
      )
    ) %>%
    ungroup()

  n_with_go <- sum(!is.na(sample_table$go_terms))
  n_with_uniprot <- sum(!is.na(sample_table$uniprot_id))

  message(
    "Sample ", sample_id,
    ": genes with GO terms = ", n_with_go,
    "; genes with UniProt IDs = ", n_with_uniprot
  )

  for (column_name in c(
    "signal_peptide_present",
    "tm_present",
    "lcr_present",
    "pfam_present",
    "pfam_terms"
  )) {
    if (!column_name %in% names(sample_table)) {
      sample_table[[column_name]] <- NA
    }
  }

  duplicate_gene_keys <- sample_table %>%
    count(sample, seq_id) %>%
    filter(n > 1L)
  if (nrow(duplicate_gene_keys) > 0L) {
    stop("Duplicated sample + seq_id keys detected for sample: ", sample_id)
  }

  compiled[[i]] <- sample_table
}

final_table <- bind_rows(compiled) %>%
  group_by(sample) %>%
  mutate(
    tAI_z = safe_zscore(tAI),
    tAI_percentile = safe_percentile(tAI),
    log_protein_length_aa = log1p(as.numeric(protein_length_aa))
  ) %>%
  ungroup() %>%
  select(
    sample,
    species,
    accession,
    domain,
    kingdom,
    phylum,
    lifestyle,
    genetic_code,
    seq_id,
    gene_id,
    protein_id,
    uniprot_id,
    protein_name,
    cds_length_nt,
    protein_length_aa,
    log_protein_length_aa,
    cds_length_multiple_of_three,
    has_terminal_stop,
    cds_qc_pass,
    metrics_available,
    trna_qc_status,
    trna_qc_reasons,
    trna_qc_pass,
    n_total_trnas,
    n_elongator_trnas,
    n_unique_anticodons,
    n_trna_amino_acids,
    any_of(c("ENC", "CAI", "FOP", "tAI", "tAI_z", "tAI_percentile", "GC", "GC3s")),
    signal_peptide_present,
    tm_present,
    lcr_present,
    pfam_present,
    pfam_terms,
    go_terms,
    everything(),
    -uniprot_id_header,
    -uniprot_id_external,
    -go_terms_header,
    -go_terms_external
  )


# Double check
n_go_annotated <- sum(!is.na(final_table$go_terms))
n_unique_go <- final_table %>%
  pull(go_terms) %>%
  str_split(";", simplify = FALSE) %>%
  unlist(use.names = FALSE) %>%
  unique() %>%
  discard(~ is.na(.x) || !nzchar(.x)) %>%
  length()

message("Rows with GO terms: ", n_go_annotated)
message("Unique GO terms: ", n_unique_go)

if (isTRUE(args$require_go_terms)) {
  if (n_go_annotated < args$min_go_annotated_genes) {
    stop(
      "GO terms are required, but only ", n_go_annotated,
      " genes have valid GO annotations. ",
      "Provide a valid --annotation-table or disable GO enrichment."
    )
  }
}


dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
write_tsv(final_table, args$output, na = "NA")

message("Saved gene table: ", args$output)
message("Rows: ", nrow(final_table))
message("Samples: ", n_distinct(final_table$sample))
message("Rows passing CDS QC: ", sum(final_table$cds_qc_pass, na.rm = TRUE))
message("Rows with finite tAI: ", sum(is.finite(final_table$tAI)))
message("Rows with GO terms: ", sum(!is.na(final_table$go_terms)))
