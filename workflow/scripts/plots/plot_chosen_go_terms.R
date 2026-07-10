#!/usr/bin/env Rscript
# Robust chosen GO term plots. The enrichment dotplot is produced only when an
# interpretable effect column can be detected. Odds-ratio columns are converted
# to log2(OR). If no effect column exists, a diagnostic file is written and a
# significance-only dotplot is produced instead.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    chosen_go_terms = "results/statistics/chosen_GOterms.tsv",
    gene_table = "results/tables/gene_features.tsv",
    output_dir = "results/plots/go_chosen_terms",
    formats = "png,pdf",
    max_terms = "48"
  )
  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--")) stop("Unexpected argument: ", args[[i]])
    if (i == length(args)) stop("Missing value for argument: ", args[[i]])
    out[[gsub("-", "_", sub("^--", "", args[[i]]))]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args()
formats <- strsplit(args$formats, ",", fixed = TRUE)[[1]] |> trimws()
max_terms <- as.integer(args$max_terms)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot_pair <- function(p, stem, width = 9, height = 7) {
  for (fmt in formats) {
    ggplot2::ggsave(
      file.path(args$output_dir, paste0(stem, ".", fmt)),
      p, width = width, height = height, units = "in", dpi = 300, bg = "white"
    )
  }
}

clean_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x, fixed = TRUE)
  stringr::str_squish(x)
}

first_existing <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) NA_character_ else hit[[1]]
}

as_num <- function(x) suppressWarnings(as.numeric(x))

if (!file.exists(args$chosen_go_terms)) stop("chosen_GOterms.tsv does not exist: ", args$chosen_go_terms)
chosen <- readr::read_tsv(args$chosen_go_terms, show_col_types = FALSE, progress = FALSE)

missing <- tibble(plot = character(), reason = character())
add_missing <- function(plot, reason) {
  missing <<- bind_rows(missing, tibble(plot = plot, reason = reason))
  message("XXXXX skipped ", plot, ": ", reason)
}

if (nrow(chosen) == 0L) {
  add_missing("all_go_plots", "chosen_GOterms.tsv has zero rows")
  readr::write_tsv(missing, file.path(args$output_dir, "XXXXX_missing_go_plots.tsv"))
  quit(save = "no", status = 0)
}

go_col <- first_existing(chosen, c("go_id", "go_term", "go", "term", "go_terms"))
name_col <- first_existing(chosen, c("go_name", "term_name", "name", "description"))
namespace_col <- first_existing(chosen, c("go_namespace", "namespace", "ontology"))
tail_col <- first_existing(chosen, c("tail", "tai_tail", "direction", "contrast", "comparison", "group"))
q_col <- first_existing(chosen, c("q_value", "q", "fdr", "padj", "adjusted_p_value", "p_adjust", "p_adj"))
p_col <- first_existing(chosen, c("p_value", "p", "pvalue"))

if (is.na(go_col)) stop("chosen_GOterms.tsv must contain a GO ID column, e.g. go_id")

log2_cols <- c("log2_odds_ratio", "log2_or", "log2OR", "log2_common_odds_ratio", "cmh_log2_or", "common_log2_odds_ratio")
log_cols <- c("log_odds_ratio", "log_or", "common_log_odds_ratio", "cmh_log_or")
or_cols <- c("odds_ratio", "OR", "common_odds_ratio", "cmh_odds_ratio", "common_or", "estimate")

effect_col <- first_existing(chosen, log2_cols)
effect_type <- "log2"
if (is.na(effect_col)) { effect_col <- first_existing(chosen, log_cols); effect_type <- "log" }
if (is.na(effect_col)) { effect_col <- first_existing(chosen, or_cols); effect_type <- "or" }

plot_df <- chosen %>%
  mutate(
    go_id_plot = as.character(.data[[go_col]]),
    go_name_plot = if (!is.na(name_col)) as.character(.data[[name_col]]) else go_id_plot,
    namespace_plot = if (!is.na(namespace_col)) clean_label(.data[[namespace_col]]) else "GO",
    tail_plot = if (!is.na(tail_col)) clean_label(.data[[tail_col]]) else "Selected GO terms",
    q_plot = if (!is.na(q_col)) as_num(.data[[q_col]]) else if (!is.na(p_col)) as_num(.data[[p_col]]) else NA_real_,
    p_plot = if (!is.na(p_col)) as_num(.data[[p_col]]) else NA_real_,
    term_label = paste0(go_id_plot, " — ", stringr::str_wrap(go_name_plot, 45))
  )

if (!is.na(effect_col)) {
  raw_effect <- as_num(chosen[[effect_col]])
  if (effect_type == "log2") plot_df$log2_effect <- raw_effect
  if (effect_type == "log") plot_df$log2_effect <- raw_effect / log(2)
  if (effect_type == "or") plot_df$log2_effect <- ifelse(raw_effect > 0, log2(raw_effect), NA_real_)
} else {
  plot_df$log2_effect <- NA_real_
}

readr::write_tsv(
  tibble(
    detected_go_column = go_col,
    detected_name_column = name_col,
    detected_namespace_column = namespace_col,
    detected_tail_column = tail_col,
    detected_q_column = q_col,
    detected_p_column = p_col,
    detected_effect_column = effect_col,
    detected_effect_type = effect_type
  ),
  file.path(args$output_dir, "chosen_go_terms_detected_columns.tsv")
)

plot_df <- plot_df %>%
  mutate(rank_score = ifelse(is.finite(q_plot) & q_plot > 0, -log10(q_plot), NA_real_)) %>%
  arrange(q_plot, desc(abs(log2_effect))) %>%
  slice_head(n = max_terms)

if (any(is.finite(plot_df$log2_effect))) {
  d <- plot_df %>% filter(is.finite(log2_effect))
  d$term_label <- factor(d$term_label, levels = rev(unique(d$term_label)))
  p <- ggplot(d, aes(x = log2_effect, y = term_label, size = rank_score, colour = tail_plot)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.85) +
    facet_grid(namespace_plot ~ ., scales = "free_y", space = "free_y") +
    labs(
      x = "log2 odds ratio / effect size",
      y = NULL,
      colour = NULL,
      size = "-log10(q/p)",
      title = "Selected GO terms from tAI-tail enrichment",
      subtitle = paste0("Effect column detected: ", effect_col, " [", effect_type, "]")
    ) +
    theme_minimal(base_size = 9) +
    theme(strip.text.y = element_text(angle = 0), legend.position = "right")
  save_plot_pair(p, "chosen_go_terms_enrichment_dotplot", width = 11, height = max(6, 0.22 * nrow(d) + 2.5))
} else {
  add_missing("chosen_go_terms_enrichment_dotplot", "no finite log2 odds/effect column found")
  d <- plot_df %>% filter(is.finite(rank_score))
  if (nrow(d) > 0L) {
    d$term_label <- factor(d$term_label, levels = rev(unique(d$term_label)))
    p <- ggplot(d, aes(x = rank_score, y = term_label, colour = tail_plot)) +
      geom_point(size = 2.3, alpha = 0.85) +
      facet_grid(namespace_plot ~ ., scales = "free_y", space = "free_y") +
      labs(
        x = "-log10(q-value or p-value)", y = NULL, colour = NULL,
        title = "Selected GO terms ranked by significance",
        subtitle = "No odds-ratio/effect column was detected, so this is not an enrichment-effect dotplot."
      ) +
      theme_minimal(base_size = 9) +
      theme(strip.text.y = element_text(angle = 0), legend.position = "right")
    save_plot_pair(p, "chosen_go_terms_significance_dotplot", width = 11, height = max(6, 0.22 * nrow(d) + 2.5))
  }
}

if (nrow(missing) == 0L) missing <- tibble(plot = "none", reason = "all requested GO plots generated")
readr::write_tsv(missing, file.path(args$output_dir, "XXXXX_missing_go_plots.tsv"))
message("Done: ", args$output_dir)
