#!/usr/bin/env Rscript
# Selects a biologically defensible set of GO terms for plotting.
#
# The script does not hand-pick GO terms from names. It starts from the actual
# genome-stratified CMH enrichment result and keeps only terms that pass explicit
# statistical and coverage filters. This makes chosen_GOterms.tsv a documented
# consequence of the analysis, not an arbitrary aesthetic selection.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

# -----------------------------------------------------------------------------
# Minimal dependency-free CLI parser. This avoids requiring r-argparse.
# -----------------------------------------------------------------------------
parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    name <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[name]] <- TRUE
      i <- i + 1L
    } else {
      out[[name]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

get_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (is.null(value)) {
    if (required) stop("Missing required argument --", name)
    return(default)
  }
  value
}

as_num <- function(x, default) {
  if (is.null(x)) return(default)
  value <- suppressWarnings(as.numeric(x))
  if (!is.finite(value)) default else value
}

as_int <- function(x, default) {
  as.integer(round(as_num(x, default)))
}

split_csv <- function(x, default = character()) {
  if (is.null(x) || !nzchar(x)) return(default)
  values <- trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE))
  values[nzchar(values)]
}

safe_read <- function(path, label) {
  if (!file.exists(path)) stop(label, " does not exist: ", path)
  read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

pick_first_col <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    stop("Cannot find ", label, " column. Tried: ", paste(candidates, collapse = ", "))
  }
  hit[[1]]
}

# -----------------------------------------------------------------------------
# Parse inputs and thresholds.
# -----------------------------------------------------------------------------
args <- parse_cli()

go_enrichment_path <- get_arg(args, "go-enrichment", required = TRUE)
gene_features_path <- get_arg(args, "gene-features", required = TRUE)
go_dictionary_path <- get_arg(args, "go-dictionary", required = TRUE)
output_path <- get_arg(args, "output", required = TRUE)
diagnostics_path <- get_arg(args, "diagnostics-output", default = sub("\\.tsv$", "_diagnostics.tsv", output_path))

q_threshold <- as_num(get_arg(args, "q-threshold", default = "0.05"), 0.05)
min_abs_log2_or <- as_num(get_arg(args, "min-abs-log2-or", default = "0.25"), 0.25)
min_informative_genomes <- as_int(get_arg(args, "min-informative-genomes", default = "5"), 5L)
min_total_genes_with_term <- as_int(get_arg(args, "min-total-genes-with-term", default = "30"), 30L)
max_terms_per_tail_namespace <- as_int(get_arg(args, "max-terms-per-tail-namespace", default = "8"), 8L)
allowed_namespaces <- split_csv(
  get_arg(args, "namespaces", default = "biological_process,molecular_function,cellular_component"),
  default = c("biological_process", "molecular_function", "cellular_component")
)

# -----------------------------------------------------------------------------
# Load and normalize enrichment result.
# -----------------------------------------------------------------------------
go <- safe_read(go_enrichment_path, "GO enrichment table")
features <- safe_read(gene_features_path, "gene_features table")
dictionary <- safe_read(go_dictionary_path, "GO dictionary")

required_go <- c("tail", "go_id", "common_odds_ratio")
missing_go <- setdiff(required_go, names(go))
if (length(missing_go) > 0L) {
  stop("GO enrichment table lacks columns: ", paste(missing_go, collapse = ", "))
}

q_col <- pick_first_col(go, c("q_value", "padj", "fdr", "adjusted_p_value"), "adjusted p-value")
p_col <- pick_first_col(go, c("p_value", "p", "raw_p_value"), "raw p-value")

if (!all(c("go_terms", "go_name", "go_namespace") %in% names(dictionary))) {
  stop("GO dictionary must contain go_terms, go_name, go_namespace")
}

# Gene-level coverage is used only as an extra sanity check/diagnostic. The
# enrichment table already contains total_genes_with_term from the CMH analysis.
if (!all(c("sample", "gene_id", "go_terms") %in% names(features))) {
  stop("gene_features must contain sample, gene_id and go_terms")
}

coverage <- features %>%
  transmute(
    sample = as.character(sample),
    gene_id = as.character(gene_id),
    go_id = str_extract_all(toupper(as.character(go_terms)), "GO:[0-9]{7}")
  ) %>%
  unnest_longer(go_id) %>%
  filter(!is.na(go_id), nzchar(go_id)) %>%
  distinct(sample, gene_id, go_id) %>%
  group_by(go_id) %>%
  summarise(
    n_genes_with_go_in_gene_features = n_distinct(paste(sample, gene_id, sep = "\t")),
    n_samples_with_go_in_gene_features = n_distinct(sample),
    .groups = "drop"
  )

dictionary_norm <- dictionary %>%
  transmute(
    go_id = as.character(go_terms),
    go_name_dictionary = as.character(go_name),
    go_namespace_dictionary = as.character(go_namespace)
  )

candidate_table <- go %>%
  mutate(
    go_id = as.character(go_id),
    tail = as.character(tail),
    q_value_selected = suppressWarnings(as.numeric(.data[[q_col]])),
    p_value_selected = suppressWarnings(as.numeric(.data[[p_col]])),
    common_odds_ratio = suppressWarnings(as.numeric(common_odds_ratio)),
    log2_common_odds_ratio = log2(common_odds_ratio),
    abs_log2_common_odds_ratio = abs(log2_common_odds_ratio),
    n_informative_genomes = suppressWarnings(as.integer(n_informative_genomes)),
    total_genes_with_term = suppressWarnings(as.integer(total_genes_with_term)),
    status = if ("status" %in% names(.)) as.character(status) else "ok"
  ) %>%
  left_join(dictionary_norm, by = "go_id") %>%
  mutate(
    go_name = coalesce(as.character(.data[[if ("go_name" %in% names(.)) "go_name" else "go_name_dictionary"]]), go_name_dictionary, go_id),
    go_namespace = coalesce(as.character(.data[[if ("go_namespace" %in% names(.)) "go_namespace" else "go_namespace_dictionary"]]), go_namespace_dictionary, "unknown")
  ) %>%
  left_join(coverage, by = "go_id") %>%
  mutate(
    n_genes_with_go_in_gene_features = coalesce(n_genes_with_go_in_gene_features, 0L),
    n_samples_with_go_in_gene_features = coalesce(n_samples_with_go_in_gene_features, 0L),
    pass_status = is.na(status) | status == "ok",
    pass_q = is.finite(q_value_selected) & q_value_selected <= q_threshold,
    pass_effect = is.finite(abs_log2_common_odds_ratio) & abs_log2_common_odds_ratio >= min_abs_log2_or,
    pass_informative_genomes = is.finite(n_informative_genomes) & n_informative_genomes >= min_informative_genomes,
    pass_total_genes = is.finite(total_genes_with_term) & total_genes_with_term >= min_total_genes_with_term,
    pass_namespace = go_namespace %in% allowed_namespaces,
    eligible = pass_status & pass_q & pass_effect & pass_informative_genomes & pass_total_genes & pass_namespace,
    enrichment_direction = case_when(
      common_odds_ratio > 1 ~ "enriched_in_tail",
      common_odds_ratio < 1 ~ "depleted_from_tail",
      TRUE ~ "neutral"
    ),
    biological_axis = case_when(
      tail == "high_tAI" & common_odds_ratio > 1 ~ "high_tAI_enriched",
      tail == "high_tAI" & common_odds_ratio < 1 ~ "high_tAI_depleted",
      tail == "low_tAI" & common_odds_ratio > 1 ~ "low_tAI_enriched",
      tail == "low_tAI" & common_odds_ratio < 1 ~ "low_tAI_depleted",
      TRUE ~ paste(tail, enrichment_direction, sep = "_")
    )
  )

chosen <- candidate_table %>%
  filter(eligible) %>%
  group_by(tail, go_namespace) %>%
  arrange(q_value_selected, desc(abs_log2_common_odds_ratio), desc(n_informative_genomes), .by_group = TRUE) %>%
  slice_head(n = max_terms_per_tail_namespace) %>%
  ungroup() %>%
  arrange(tail, go_namespace, q_value_selected, desc(abs_log2_common_odds_ratio)) %>%
  transmute(
    go_id,
    go_name,
    go_namespace,
    tail,
    biological_axis,
    enrichment_direction,
    common_odds_ratio,
    log2_common_odds_ratio,
    p_value = p_value_selected,
    q_value = q_value_selected,
    n_genomes_with_term,
    n_informative_genomes,
    total_genes_with_term,
    n_genes_with_go_in_gene_features,
    n_samples_with_go_in_gene_features,
    selection_reason = paste0(
      "CMH tail enrichment; q<=", q_threshold,
      "; abs(log2OR)>=", min_abs_log2_or,
      "; informative_genomes>=", min_informative_genomes,
      "; total_genes_with_term>=", min_total_genes_with_term
    )
  )

# If no terms pass, write an explicit placeholder rather than silently selecting
# weak terms. This follows the project rule: do not generate biological plots
# from uncertain information.
if (nrow(chosen) == 0L) {
  chosen <- tibble(
    go_id = "XXXXX",
    go_name = "XXXXX",
    go_namespace = "XXXXX",
    tail = "XXXXX",
    biological_axis = "XXXXX",
    enrichment_direction = "XXXXX",
    common_odds_ratio = NA_real_,
    log2_common_odds_ratio = NA_real_,
    p_value = NA_real_,
    q_value = NA_real_,
    n_genomes_with_term = NA_integer_,
    n_informative_genomes = NA_integer_,
    total_genes_with_term = NA_integer_,
    n_genes_with_go_in_gene_features = NA_integer_,
    n_samples_with_go_in_gene_features = NA_integer_,
    selection_reason = "XXXXX: no GO term passed the configured biological/statistical filters"
  )
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(diagnostics_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(chosen, output_path, na = "NA")
write_tsv(candidate_table, diagnostics_path, na = "NA")

message("Saved chosen GO terms: ", output_path)
message("Saved GO term diagnostics: ", diagnostics_path)
message("Chosen rows: ", nrow(chosen))
