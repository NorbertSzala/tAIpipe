#!/usr/bin/env Rscript
# Tests PFAM-domain enrichment in the within-genome top and bottom 1%/10% tAI
# tails. Genome-stratified CMH tests prevent large genomes from dominating the
# effect estimate; the output also retains counts and mean raw tAI for reporting.

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

parser <- ArgumentParser(description = "PFAM enrichment in within-genome tAI tails")
parser$add_argument("--gene-features", required = TRUE)
parser$add_argument("--pfam-description-table", default = "")
parser$add_argument("--tail-fractions", default = "0.10,0.01")
parser$add_argument("--min-genes-per-sample", type = "integer", default = 100L)
parser$add_argument("--min-tail-size", type = "integer", default = 20L)
parser$add_argument("--min-informative-genomes", type = "integer", default = 3L)
parser$add_argument("--min-total-genes-with-term", type = "integer", default = 10L)
parser$add_argument("--fdr-method", default = "BH")
parser$add_argument("--output", required = TRUE)
parser$add_argument("--qc-output", required = TRUE)
args <- parser$parse_args()

tail_fractions <- as.numeric(strsplit(args$tail_fractions, ",", fixed = TRUE)[[1]])
if (
  length(tail_fractions) == 0L || any(!is.finite(tail_fractions)) ||
  any(tail_fractions <= 0 | tail_fractions >= 0.5)
) {
  stop("--tail-fractions must contain comma-separated values in (0, 0.5)")
}
if (args$min_genes_per_sample < 2L || args$min_tail_size < 1L) {
  stop("Minimum sample and tail sizes must be positive")
}

normalize_pfam_id <- function(x) {
  stringr::str_extract(toupper(as.character(x)), "PF[0-9]{5}")
}

extract_pfam_ids <- function(x) {
  ids <- stringr::str_extract_all(toupper(as.character(x)), "PF[0-9]{5}")[[1]]
  unique(ids[!is.na(ids) & nzchar(ids)])
}

read_pfam_descriptions <- function(path) {
  empty <- tibble(PFAM = character(), description = character())
  if (!nzchar(path) || !file.exists(path)) return(empty)
  tab <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  id_col <- intersect(c("PFAM", "pfam_id", "accession", "pfam_acc", "AC", "ACC"), names(tab))
  desc_col <- intersect(c("description", "pfam_description", "pfam_name", "name", "function", "DE", "DESC", "NAME"), names(tab))
  if (length(id_col) == 0L || length(desc_col) == 0L) {
    warning("PFAM description table lacks recognized identifier/description columns: ", path)
    return(empty)
  }
  tab %>%
    transmute(
      PFAM = normalize_pfam_id(.data[[id_col[[1]]]]),
      description = stringr::str_squish(stringr::str_replace_all(
        as.character(.data[[desc_col[[1]]]]), "[\r\n\t]+", " "
      ))
    ) %>%
    filter(!is.na(PFAM), nzchar(description)) %>%
    distinct(PFAM, .keep_all = TRUE)
}

safe_cmh <- function(strata) {
  k <- nrow(strata)
  arr <- array(0, dim = c(2L, 2L, k))
  for (i in seq_len(k)) {
    arr[, , i] <- matrix(
      c(
        strata$tail_pfam[[i]], strata$tail_without_pfam[[i]],
        strata$background_pfam[[i]], strata$background_without_pfam[[i]]
      ),
      nrow = 2L,
      byrow = TRUE
    )
  }
  fit <- tryCatch(
    stats::mantelhaen.test(arr, correct = FALSE, exact = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble(
      common_odds_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
      statistic = NA_real_, p_value = NA_real_, status = paste0("cmh_error: ", conditionMessage(fit))
    ))
  }
  tibble(
    common_odds_ratio = unname(fit$estimate),
    conf_low = unname(fit$conf.int[[1]]),
    conf_high = unname(fit$conf.int[[2]]),
    statistic = unname(fit$statistic),
    p_value = fit$p.value,
    status = "ok"
  )
}

raw <- readr::read_tsv(args$gene_features, show_col_types = FALSE, progress = FALSE)
required <- c("sample", "tAI", "pfam_terms")
missing <- setdiff(required, names(raw))
if (length(missing) > 0L) stop("gene_features.tsv lacks columns: ", paste(missing, collapse = ", "))
if (!"gene_id" %in% names(raw)) {
  raw$gene_id <- if ("seq_id" %in% names(raw)) as.character(raw$seq_id) else paste0("gene_", seq_len(nrow(raw)))
}

genes <- raw %>%
  transmute(
    sample = as.character(sample),
    gene_id = as.character(gene_id),
    tAI = suppressWarnings(as.numeric(tAI)),
    pfam_terms = as.character(pfam_terms)
  ) %>%
  filter(!is.na(sample), nzchar(sample), !is.na(gene_id), nzchar(gene_id), is.finite(tAI)) %>%
  distinct(sample, gene_id, .keep_all = TRUE)
if (nrow(genes) == 0L) stop("No genes with finite tAI are available")

pfam_long <- genes %>%
  mutate(PFAM = purrr::map(pfam_terms, extract_pfam_ids)) %>%
  select(sample, gene_id, tAI, PFAM) %>%
  tidyr::unnest_longer(PFAM) %>%
  filter(!is.na(PFAM), nzchar(PFAM)) %>%
  distinct(sample, gene_id, PFAM, .keep_all = TRUE)
if (nrow(pfam_long) == 0L) stop("No PFAM identifiers were extracted from pfam_terms")

select_membership <- function(fraction) {
  genes %>%
    group_by(sample) %>%
    arrange(tAI, gene_id, .by_group = TRUE) %>%
    mutate(
      n_genes = n(),
      tail_size = floor(n_genes * fraction),
      eligible = n_genes >= args$min_genes_per_sample & tail_size >= args$min_tail_size,
      rank_low = row_number(),
      rank_high = n_genes - row_number() + 1L,
      low_tAI = eligible & rank_low <= tail_size,
      high_tAI = eligible & rank_high <= tail_size,
      tail_fraction = fraction
    ) %>%
    ungroup()
}

memberships <- purrr::map(tail_fractions, select_membership)
qc <- purrr::map_dfr(memberships, function(x) {
  x %>%
    group_by(sample, tail_fraction) %>%
    summarise(
      n_genes_finite_tAI = first(n_genes),
      n_tail_each_direction = first(tail_size),
      eligible = first(eligible),
      .groups = "drop"
    )
})
dir.create(dirname(args$qc_output), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(qc, args$qc_output, na = "NA")

analyse_direction <- function(membership, direction) {
  fraction <- unique(membership$tail_fraction)[[1]]
  membership <- membership %>% filter(eligible) %>%
    transmute(sample, gene_id, tAI, is_tail = .data[[direction]])
  if (nrow(membership) == 0L) return(tibble())

  totals <- membership %>%
    group_by(sample) %>%
    summarise(tail_total = sum(is_tail), background_total = sum(!is_tail), .groups = "drop")

  counts <- pfam_long %>%
    inner_join(membership %>% select(sample, gene_id, is_tail), by = c("sample", "gene_id")) %>%
    group_by(sample, PFAM) %>%
    summarise(
      tail_pfam = n_distinct(gene_id[is_tail]),
      background_pfam = n_distinct(gene_id[!is_tail]),
      .groups = "drop"
    ) %>%
    left_join(totals, by = "sample") %>%
    mutate(
      tail_without_pfam = tail_total - tail_pfam,
      background_without_pfam = background_total - background_pfam,
      genes_with_term = tail_pfam + background_pfam,
      genes_without_term = tail_without_pfam + background_without_pfam,
      informative = tail_total > 0L & background_total > 0L & genes_with_term > 0L & genes_without_term > 0L,
      stratum_log_odds_ratio = log(
        ((tail_pfam + 0.5) * (background_without_pfam + 0.5)) /
          ((tail_without_pfam + 0.5) * (background_pfam + 0.5))
      )
    )

  metadata <- counts %>%
    group_by(PFAM) %>%
    summarise(
      n_count = sum(tail_pfam),
      background_count = sum(background_pfam),
      n_genomes_with_term = n_distinct(sample),
      n_genomes_with_term_in_tail = n_distinct(sample[tail_pfam > 0]),
      n_informative_genomes = n_distinct(sample[informative]),
      total_genes_with_term = sum(genes_with_term),
      n_genomes_enriched = sum(informative & stratum_log_odds_ratio > 0),
      n_genomes_depleted = sum(informative & stratum_log_odds_ratio < 0),
      .groups = "drop"
    ) %>%
    filter(
      n_informative_genomes >= args$min_informative_genomes,
      total_genes_with_term >= args$min_total_genes_with_term,
      n_count > 0L
    )
  if (nrow(metadata) == 0L) return(tibble())

  means <- pfam_long %>%
    inner_join(membership, by = c("sample", "gene_id", "tAI")) %>%
    filter(is_tail) %>%
    group_by(PFAM) %>%
    summarise(mean_tAI = mean(tAI, na.rm = TRUE), median_tAI = median(tAI, na.rm = TRUE), .groups = "drop")

  metadata %>%
    split(.$PFAM) %>%
    purrr::map_dfr(function(meta) {
      id <- meta$PFAM[[1]]
      strata <- counts %>% filter(PFAM == id, informative) %>% arrange(sample)
      bind_cols(meta, safe_cmh(strata))
    }) %>%
    left_join(means, by = "PFAM") %>%
    mutate(
      tail_fraction = fraction,
      tail_percent = 100 * fraction,
      direction = direction,
      label = if_else(direction == "high_tAI", "highest", "lowest"),
      tail_definition = paste0(if_else(direction == "high_tAI", "Top ", "Bottom "), format(100 * fraction, trim = TRUE), "% tAI within each genome")
    )
}

results <- purrr::map_dfr(memberships, function(membership) {
  bind_rows(
    analyse_direction(membership, "high_tAI"),
    analyse_direction(membership, "low_tAI")
  )
})

descriptions <- read_pfam_descriptions(args$pfam_description_table)
if (nrow(results) == 0L) {
  results <- tibble(
    PFAM = character(), description = character(), n_count = integer(),
    enrichment_overall = double(), mean_tAI = double(), label = character(),
    source_df = character(), additional = character(), tail_fraction = double(),
    tail_percent = double(), tail_definition = character(), direction = character(),
    q_value = double()
  )
} else {
  results <- results %>%
    group_by(tail_fraction, direction) %>%
    mutate(q_value = p.adjust(p_value, method = args$fdr_method)) %>%
    ungroup() %>%
    left_join(descriptions, by = "PFAM") %>%
    mutate(
      description = dplyr::coalesce(description, "Description not supplied"),
      enrichment_overall = common_odds_ratio,
      log2_enrichment = if_else(common_odds_ratio > 0, log2(common_odds_ratio), -Inf),
      source_df = "gene_features.tsv; PFAM tail vs within-genome remainder",
      additional = NA_character_
    ) %>%
    select(
      PFAM, description, n_count, enrichment_overall, mean_tAI, label,
      source_df, additional, tail_fraction, tail_percent, tail_definition,
      direction, background_count, median_tAI, n_genomes_with_term,
      n_genomes_with_term_in_tail, n_informative_genomes, total_genes_with_term,
      n_genomes_enriched, n_genomes_depleted, common_odds_ratio, log2_enrichment,
      conf_low, conf_high, statistic, p_value, q_value, status
    ) %>%
    arrange(tail_fraction, direction, q_value, desc(abs(log2_enrichment)), desc(n_count))
}

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(results, args$output, na = "NA")
message("Saved PFAM tAI-tail enrichment: ", args$output)
message("Rows: ", nrow(results))

