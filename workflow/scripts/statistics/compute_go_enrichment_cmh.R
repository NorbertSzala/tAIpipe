#!/usr/bin/env Rscript
# Performs genome-stratified GO enrichment analysis for the high- and low-tAI
# tails. Tails are selected independently within each genome from the
# GO-annotated gene universe, which avoids treating annotation absence as true
# GO-term absence.

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

parser <- ArgumentParser(
  description = paste(
    "Run GO enrichment for high- and low-tAI tails using a",
    "Cochran-Mantel-Haenszel test stratified by genome."
  )
)
parser$add_argument("--gene-features", required = TRUE)
parser$add_argument("--go-dictionary", required = TRUE)
parser$add_argument("--output", required = TRUE)
parser$add_argument("--tail-fraction", type = "double", default = 0.10)
parser$add_argument(
  "--max-tail-genes",
  type = "integer",
  default = 0L,
  help = "Zero disables the absolute cap."
)
parser$add_argument("--min-genomes-with-term", type = "integer", default = 5L)
parser$add_argument("--min-informative-genomes", type = "integer", default = 3L)
parser$add_argument("--min-total-genes-with-term", type = "integer", default = 10L)
parser$add_argument("--min-genes-with-go-total", type = "integer", default = 100L)
parser$add_argument("--min-genes-with-go-per-sample", type = "integer", default = 10L)
parser$add_argument("--fdr-method", default = "BH")
args <- parser$parse_args()

if (!is.finite(args$tail_fraction) || args$tail_fraction <= 0 || args$tail_fraction >= 0.5) {
  stop("--tail-fraction must be greater than 0 and lower than 0.5")
}
if (args$max_tail_genes < 0L) {
  stop("--max-tail-genes must be zero or a positive integer")
}
if (args$min_genes_with_go_total < 1L) {
  stop("--min-genes-with-go-total must be a positive integer")
}
if (args$min_genes_with_go_per_sample < 1L) {
  stop("--min-genes-with-go-per-sample must be a positive integer")
}

extract_go_ids <- function(x) {
  ids <- str_extract_all(toupper(as.character(x)), "GO:[0-9]{7}")[[1]]
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  ids
}

required_columns <- c("sample", "tAI", "go_terms")
raw_genes <- read_tsv(args$gene_features, show_col_types = FALSE)
missing_columns <- setdiff(required_columns, names(raw_genes))
if (length(missing_columns) > 0L) {
  stop("Gene feature table lacks columns: ", paste(missing_columns, collapse = ", "))
}

if (!"gene_id" %in% names(raw_genes)) {
  if ("seq_id" %in% names(raw_genes)) {
    raw_genes <- raw_genes %>% mutate(gene_id = as.character(seq_id))
  } else {
    raw_genes <- raw_genes %>% mutate(gene_id = paste0("gene_", row_number()))
  }
}

genes <- raw_genes %>%
  transmute(
    sample = as.character(sample),
    gene_id = as.character(gene_id),
    tAI = as.numeric(tAI),
    go_terms = as.character(go_terms)
  ) %>%
  filter(
    !is.na(sample), sample != "",
    !is.na(gene_id), gene_id != "",
    is.finite(tAI)
  ) %>%
  distinct(sample, gene_id, .keep_all = TRUE) %>%
  mutate(
    go_id_list = map(go_terms, extract_go_ids),
    has_go = lengths(go_id_list) > 0L
  )

if (nrow(genes) == 0L) {
  stop("No genes with finite tAI are available.")
}

n_annotated_total <- sum(genes$has_go)
if (n_annotated_total < args$min_genes_with_go_total) {
  stop(
    "GO enrichment requires at least ", args$min_genes_with_go_total,
    " genes with valid GO terms after tAI filtering, but found ",
    n_annotated_total, "."
  )
}

sample_annotation_qc <- genes %>%
  group_by(sample) %>%
  summarise(
    n_genes_finite_tai = n(),
    n_genes_with_go = sum(has_go),
    go_annotated_fraction = n_genes_with_go / n_genes_finite_tai,
    .groups = "drop"
  ) %>%
  mutate(
    sample_used_for_go_enrichment = n_genes_with_go >= args$min_genes_with_go_per_sample
  )

eligible_samples <- sample_annotation_qc %>%
  filter(sample_used_for_go_enrichment) %>%
  pull(sample)

if (length(eligible_samples) == 0L) {
  stop(
    "No sample has at least ", args$min_genes_with_go_per_sample,
    " genes with valid GO terms."
  )
}

ranked <- genes %>%
  filter(sample %in% eligible_samples, has_go) %>%
  select(sample, gene_id, tAI, go_id_list) %>%
  group_by(sample) %>%
  arrange(desc(tAI), gene_id, .by_group = TRUE) %>%
  mutate(rank_high = row_number()) %>%
  arrange(tAI, gene_id, .by_group = TRUE) %>%
  mutate(
    rank_low = row_number(),
    n_eligible_genes = n(),
    n_tail_fraction = ceiling(args$tail_fraction * n_eligible_genes),
    n_tail_uncapped = if (args$max_tail_genes > 0L) {
      pmin(n_tail_fraction, args$max_tail_genes)
    } else {
      n_tail_fraction
    },
    # Prevent high and low tails from overlapping in very small genomes.
    n_tail = pmin(n_tail_uncapped, floor(n_eligible_genes / 2)),
    high_tail = rank_high <= n_tail,
    low_tail = rank_low <= n_tail
  ) %>%
  ungroup()

if (nrow(ranked) == 0L) {
  stop("No GO-annotated genes remained after sample-level GO coverage filtering.")
}

small_tail_samples <- ranked %>%
  group_by(sample) %>%
  summarise(n_tail = min(n_tail), n_eligible_genes = n(), .groups = "drop") %>%
  filter(n_tail < 1L)
if (nrow(small_tail_samples) > 0L) {
  stop("At least one eligible sample has no selectable tail genes after filtering.")
}

go_long <- ranked %>%
  unnest_longer(go_id_list, values_to = "go_id") %>%
  filter(!is.na(go_id), go_id != "") %>%
  distinct(sample, gene_id, go_id, high_tail, low_tail)

if (nrow(go_long) == 0L) {
  stop("No valid GO identifiers were extracted from the go_terms column.")
}

safe_cmh <- function(strata) {
  k <- nrow(strata)
  array_2x2xk <- array(0, dim = c(2L, 2L, k))

  for (i in seq_len(k)) {
    array_2x2xk[, , i] <- matrix(
      c(
        strata$tail_go[[i]],
        strata$tail_without_go[[i]],
        strata$background_go[[i]],
        strata$background_without_go[[i]]
      ),
      nrow = 2L,
      byrow = TRUE
    )
  }

  result <- tryCatch(
    mantelhaen.test(array_2x2xk, correct = FALSE, exact = FALSE),
    error = function(error) error
  )

  if (inherits(result, "error")) {
    return(tibble(
      common_odds_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      status = paste0("cmh_error: ", conditionMessage(result))
    ))
  }

  tibble(
    common_odds_ratio = unname(result$estimate),
    conf_low = unname(result$conf.int[[1]]),
    conf_high = unname(result$conf.int[[2]]),
    statistic = unname(result$statistic),
    p_value = result$p.value,
    status = "ok"
  )
}

analyse_tail <- function(tail_column, tail_name) {
  membership <- ranked %>%
    transmute(
      sample,
      gene_id,
      is_tail = .data[[tail_column]]
    )

  totals <- membership %>%
    group_by(sample) %>%
    summarise(
      tail_total = sum(is_tail),
      background_total = sum(!is_tail),
      .groups = "drop"
    )

  counts <- go_long %>%
    transmute(
      sample,
      gene_id,
      go_id,
      is_tail = .data[[tail_column]]
    ) %>%
    distinct() %>%
    group_by(sample, go_id) %>%
    summarise(
      tail_go = n_distinct(gene_id[is_tail]),
      background_go = n_distinct(gene_id[!is_tail]),
      .groups = "drop"
    ) %>%
    left_join(totals, by = "sample") %>%
    mutate(
      tail_without_go = tail_total - tail_go,
      background_without_go = background_total - background_go,
      genes_with_term = tail_go + background_go,
      genes_without_term = tail_without_go + background_without_go,
      informative = (
        tail_total > 0L &
          background_total > 0L &
          genes_with_term > 0L &
          genes_without_term > 0L
      ),
      stratum_log_odds_ratio = log(
        ((tail_go + 0.5) * (background_without_go + 0.5)) /
          ((tail_without_go + 0.5) * (background_go + 0.5))
      )
    )

  term_metadata <- counts %>%
    group_by(go_id) %>%
    summarise(
      n_genomes_with_term = n_distinct(sample),
      n_informative_genomes = n_distinct(sample[informative]),
      total_genes_with_term = sum(genes_with_term),
      n_genomes_enriched = sum(informative & stratum_log_odds_ratio > 0),
      n_genomes_depleted = sum(informative & stratum_log_odds_ratio < 0),
      n_genomes_neutral = sum(informative & stratum_log_odds_ratio == 0),
      .groups = "drop"
    ) %>%
    filter(
      n_genomes_with_term >= args$min_genomes_with_term,
      n_informative_genomes >= args$min_informative_genomes,
      total_genes_with_term >= args$min_total_genes_with_term
    )

  if (nrow(term_metadata) == 0L) {
    return(tibble())
  }

  term_metadata %>%
    split(.$go_id) %>%
    map_dfr(function(meta_row) {
      term_id <- meta_row$go_id[[1]]
      strata <- counts %>%
        filter(go_id == term_id, informative) %>%
        arrange(sample)

      bind_cols(
        tibble(
          tail = tail_name,
          go_id = term_id
        ),
        meta_row %>% select(-go_id),
        safe_cmh(strata)
      )
    })
}

empty_results <- tibble(
  tail = character(),
  go_id = character(),
  go_name = character(),
  go_namespace = character(),
  n_genomes_with_term = integer(),
  n_informative_genomes = integer(),
  total_genes_with_term = integer(),
  n_genomes_enriched = integer(),
  n_genomes_depleted = integer(),
  n_genomes_neutral = integer(),
  common_odds_ratio = double(),
  conf_low = double(),
  conf_high = double(),
  statistic = double(),
  p_value = double(),
  status = character(),
  q_value = double(),
  go_universe = character(),
  n_annotated_genes_total = integer(),
  n_samples_in_universe = integer()
)

results <- bind_rows(
  analyse_tail("high_tail", "high_tAI"),
  analyse_tail("low_tail", "low_tAI")
)

if (nrow(results) == 0L) {
  warning("No GO terms passed the configured minimum filters.")
  results <- empty_results %>% select(-go_name, -go_namespace)
} else {
  results <- results %>%
    group_by(tail) %>%
    mutate(q_value = p.adjust(p_value, method = args$fdr_method)) %>%
    ungroup() %>%
    mutate(
      go_universe = "go_annotated_genes_with_finite_tai",
      n_annotated_genes_total = nrow(ranked),
      n_samples_in_universe = n_distinct(ranked$sample)
    )
}

go_dictionary <- read_tsv(args$go_dictionary, show_col_types = FALSE)
if (all(c("go_terms", "go_name", "go_namespace") %in% names(go_dictionary))) {
  dictionary <- go_dictionary %>%
    transmute(
      go_id = toupper(as.character(go_terms)),
      go_name = as.character(go_name),
      go_namespace = as.character(go_namespace)
    ) %>%
    distinct(go_id, .keep_all = TRUE)

  results <- results %>%
    left_join(dictionary, by = "go_id") %>%
    relocate(go_name, go_namespace, .after = go_id)

  missing_dictionary_terms <- results %>%
    filter(is.na(go_name) | is.na(go_namespace)) %>%
    distinct(go_id) %>%
    pull(go_id)

  if (length(missing_dictionary_terms) > 0L) {
    warning(
      "GO dictionary lacks names/namespaces for ", length(missing_dictionary_terms),
      " enriched terms. Example IDs: ",
      paste(head(missing_dictionary_terms, 10), collapse = ", ")
    )
  }
} else {
  results <- results %>%
    mutate(go_name = NA_character_, go_namespace = NA_character_) %>%
    relocate(go_name, go_namespace, .after = go_id)
}

results <- results %>%
  arrange(tail, q_value, p_value, desc(abs(log(common_odds_ratio))))

qc_path <- sub("\\.tsv$", "_annotation_qc.tsv", args$output)
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
write_tsv(results, args$output, na = "NA")
write_tsv(sample_annotation_qc, qc_path, na = "NA")

message("Saved CMH GO enrichment table: ", args$output)
message("Saved GO annotation QC table: ", qc_path)
message("Rows: ", nrow(results))
message("GO universe: GO-annotated genes with finite tAI; genes = ", nrow(ranked), "; samples = ", n_distinct(ranked$sample))
