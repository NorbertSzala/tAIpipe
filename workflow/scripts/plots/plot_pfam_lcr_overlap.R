#!/usr/bin/env Rscript
# PFAM-LCR overlap plots.
# tAI_z is the default because it removes genome-specific shifts in tAI scale.
# Raw tAI plots are also created when tAI is available, but should be interpreted
# as absolute codon-adaptation distributions and are more confounded by genome.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

source("workflow/scripts/lib/plot_style_helpers.R")

parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    gene_table = "results/tables/gene_features.tsv",
    output_dir = "results/plots/pfam_lcr_overlap",
    pfam_description_table = "",
    top_n = "25",
    formats = "png,pdf",
    max_plot_rows = "160000",
    seed = "1"
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
args <- if (exists("snakemake")) {
  list(
    gene_table = snakemake@input[["gene_table"]] %||% snakemake@input[["gene_features"]],
    output_dir = snakemake@params[["output_dir"]] %||% "results/plots/pfam_lcr_overlap",
    pfam_description_table = snakemake@params[["pfam_description_table"]] %||% "",
    top_n = as.character(snakemake@params[["top_n"]] %||% 25),
    formats = paste(snakemake@params[["formats"]] %||% c("png", "pdf"), collapse = ","),
    max_plot_rows = as.character(snakemake@params[["max_plot_rows"]] %||% 160000),
    seed = as.character(snakemake@params[["seed"]] %||% 1)
  )
} else parse_cli()

output_dir <- args$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
formats <- strsplit(args$formats, ",", fixed = TRUE)[[1]] |> trimws()
top_n <- as.integer(args$top_n)
max_plot_rows <- as.integer(args$max_plot_rows)
set.seed(as.integer(args$seed))

sample_rows <- function(df, max_n) {
  # dplyr::slice_sample(n = min(n(), max_n)) is invalid because n() is not
  # allowed inside the n argument. Compute the target size outside dplyr.
  if (!is.finite(max_n) || max_n <= 0L) return(df)
  k <- min(nrow(df), max_n)
  if (k <= 0L) return(df[0, , drop = FALSE])
  if (nrow(df) <= k) return(df)
  dplyr::slice_sample(df, n = k, replace = FALSE)
}

write_missing <- function(text) {
  writeLines(text, file.path(output_dir, "XXXXX_PFAM_LCR_PLOTS_NOT_CREATED.txt"))
  message(text)
}

if (!file.exists(args$gene_table)) stop("Missing gene_features table: ", args$gene_table)
genes <- readr::read_tsv(args$gene_table, show_col_types = FALSE, progress = FALSE)

if (!"pfam_lcr_overlap_present" %in% names(genes) && "pfam_lcr_present" %in% names(genes)) {
  genes$pfam_lcr_overlap_present <- genes$pfam_lcr_present
}
if (!"pfam_lcr_overlap_terms" %in% names(genes) && "pfam_lcr_terms" %in% names(genes)) {
  genes$pfam_lcr_overlap_terms <- genes$pfam_lcr_terms
}

required <- c("sample", "tAI_z", "pfam_lcr_overlap_present")
missing <- setdiff(required, names(genes))
if (length(missing) > 0L) {
  write_missing(paste("XXXXX: missing required columns:", paste(missing, collapse = ", ")))
  quit(save = "no", status = 0)
}

if (!"gene_id" %in% names(genes)) genes$gene_id <- if ("seq_id" %in% names(genes)) genes$seq_id else seq_len(nrow(genes))
if (!"phylum" %in% names(genes)) genes$phylum <- "Unknown"
if (!"lifestyle" %in% names(genes)) genes$lifestyle <- "Unknown"
if (!"pfam_present" %in% names(genes)) genes$pfam_present <- NA
if (!"lcr_present" %in% names(genes)) genes$lcr_present <- NA
if (!"pfam_terms" %in% names(genes)) genes$pfam_terms <- NA_character_
if (!"pfam_lcr_overlap_terms" %in% names(genes)) genes$pfam_lcr_overlap_terms <- NA_character_
if (!"lcr_total_length" %in% names(genes)) genes$lcr_total_length <- NA_real_
if (!"tAI" %in% names(genes)) genes$tAI <- NA_real_

plot_data <- genes %>%
  mutate(
    tAI_z = suppressWarnings(as.numeric(tAI_z)),
    tAI = suppressWarnings(as.numeric(tAI)),
    pfam_lcr_overlap_present = coerce_binary_plot(pfam_lcr_overlap_present),
    pfam_present_bool = ifelse(!is.na(coerce_binary_plot(pfam_present)), coerce_binary_plot(pfam_present), !is.na(pfam_terms)),
    lcr_present_bool = coerce_binary_plot(lcr_present),
    lcr_total_length = suppressWarnings(as.numeric(lcr_total_length)),
    phylum = clean_plot_label(phylum),
    lifestyle = clean_plot_label(lifestyle),
    overlap_label = factor(
      if_else(pfam_lcr_overlap_present, "PFAM-LCR overlap detected", "No PFAM-LCR overlap recorded"),
      levels = c("No PFAM-LCR overlap recorded", "PFAM-LCR overlap detected")
    ),
    structural_category = case_when(
      pfam_lcr_overlap_present ~ "PFAM-LCR overlap detected",
      pfam_present_bool & lcr_present_bool ~ "PFAM and LCR, no overlap recorded",
      pfam_present_bool & !lcr_present_bool ~ "PFAM only",
      !pfam_present_bool & lcr_present_bool ~ "LCR only",
      TRUE ~ "Neither PFAM nor LCR"
    )
  )

if (sum(plot_data$pfam_lcr_overlap_present, na.rm = TRUE) < 10L) {
  write_missing("XXXXX: fewer than 10 genes have PFAM-LCR overlap; plot would not be stable.")
  quit(save = "no", status = 0)
}

plot_metric_overlap <- function(metric_col, y_label, stem) {
  if (!metric_col %in% names(plot_data)) return(NULL)
  df <- plot_data %>% filter(is.finite(.data[[metric_col]]), !is.na(overlap_label))
  if (nrow(df) < 20L || n_distinct(df$overlap_label) < 2L) return(NULL)
  per_genome <- df %>% group_by(sample, overlap_label) %>% summarise(v = median(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  p <- safe_wilcox_p(per_genome$v, per_genome$overlap_label)
  p_plot <- sample_rows(df, max_plot_rows) %>%
    ggplot(aes(x = overlap_label, y = .data[[metric_col]], fill = overlap_label)) +
    geom_violin(trim = FALSE, alpha = 0.70) +
    geom_boxplot(width = 0.12, outlier.alpha = 0.18, alpha = 0.85) +
    annotate_top_right(sig_label(p)) +
    labs(
      x = NULL, y = y_label,
      title = paste(y_label, "by PFAM-LCR overlap status"),
      subtitle = "p-value is Wilcoxon test on per-genome medians; violins use sampled genes for rendering."
    ) +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  save_plot_pair(p_plot, stem, output_dir, 7.8, 5.8, formats)
}

plot_metric_overlap("tAI_z", "within-genome standardized tAI (tAI_z)", "pfam_lcr_overlap_tai_z_violin")
plot_metric_overlap("tAI", "raw tAI", "pfam_lcr_overlap_raw_tai_violin")

# Structural categories: use violin, not boxplot only. Wrap long labels.
df_cat <- plot_data %>% filter(is.finite(tAI_z), !is.na(structural_category))
p_cat <- safe_kruskal_p(df_cat$tAI_z, df_cat$structural_category)
cat_plot <- df_cat %>%
  mutate(structural_category = stringr::str_wrap(structural_category, width = 22)) %>%
  sample_rows(max_plot_rows) %>%
  ggplot(aes(x = structural_category, y = tAI_z, fill = structural_category)) +
  geom_violin(trim = FALSE, alpha = 0.68) +
  geom_boxplot(width = 0.11, outlier.alpha = 0.18, alpha = 0.90) +
  annotate_top_right(sig_label(p_cat)) +
  labs(
    x = NULL, y = "within-genome standardized tAI (tAI_z)",
    title = "tAI across PFAM/LCR structural annotation categories",
    subtitle = "Long category labels are wrapped. This is descriptive; overlap does not prove functional disruption."
  ) +
  theme_minimal(base_size = 11) + theme(legend.position = "none")
save_plot_pair(cat_plot, "pfam_lcr_structural_categories_tai_z_violin", output_dir, 10.5, 6.2, formats)
# Backward-compatible filename requested in existing targets.
save_plot_pair(cat_plot, "pfam_lcr_structural_categories_tai_z_boxplot", output_dir, 10.5, 6.2, formats)

# Long PFAM table.
pfam_long <- plot_data %>%
  filter(pfam_lcr_overlap_present, !is.na(pfam_lcr_overlap_terms), nzchar(as.character(pfam_lcr_overlap_terms))) %>%
  transmute(sample, gene_id, phylum, lifestyle, tAI_z, pfam_id = str_split(as.character(pfam_lcr_overlap_terms), ";|\\|")) %>%
  tidyr::unnest_longer(pfam_id) %>%
  mutate(pfam_id = str_trim(pfam_id)) %>%
  filter(!is.na(pfam_id), nzchar(pfam_id), pfam_id != "NA")

pfam_desc <- tibble(pfam_id = character(), pfam_description = character())
if (nzchar(args$pfam_description_table) && file.exists(args$pfam_description_table)) {
  pfam_desc <- readr::read_tsv(args$pfam_description_table, show_col_types = FALSE)
  if (!all(c("pfam_id", "pfam_description") %in% names(pfam_desc))) {
    stop("PFAM description table must contain pfam_id and pfam_description columns")
  }
}

if (nrow(pfam_long) > 0L) {
  pfam_summary <- pfam_long %>%
    group_by(pfam_id) %>%
    summarise(
      n_genes = n_distinct(paste(sample, gene_id, sep = "\t")),
      n_samples = n_distinct(sample),
      median_tAI_z = median(tAI_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(pfam_desc, by = "pfam_id") %>%
    mutate(pfam_description = dplyr::coalesce(pfam_description, "Description not supplied")) %>%
    arrange(desc(n_genes), desc(n_samples))
  readr::write_tsv(pfam_summary, file.path(output_dir, "pfam_lcr_overlap_domain_summary.tsv"))
  top_ids <- head(pfam_summary$pfam_id, top_n)

  base_counts <- pfam_long %>% filter(pfam_id %in% top_ids)
  pfam_levels <- pfam_summary %>% filter(pfam_id %in% top_ids) %>% arrange(n_genes) %>% pull(pfam_id)

  total_plot <- pfam_summary %>% filter(pfam_id %in% top_ids) %>%
    mutate(pfam_id = factor(pfam_id, levels = pfam_levels), label = paste0(pfam_id, "\n", stringr::str_wrap(pfam_description, 45))) %>%
    ggplot(aes(x = pfam_id, y = n_genes)) + geom_col() + coord_flip() +
    scale_x_discrete(labels = setNames(paste0(pfam_summary$pfam_id, "\n", stringr::str_wrap(pfam_summary$pfam_description, 45)), pfam_summary$pfam_id)) +
    labs(x = NULL, y = "Genes with PFAM-LCR overlap", title = "Most frequent PFAM domains overlapping LCRs") +
    theme_minimal(base_size = 10)
  save_plot_pair(total_plot, "top_pfam_domains_overlapping_lcr", output_dir, 9.5, 7.2, formats)

  stacked_plot <- function(group_col, stem, title) {
    d <- base_counts %>%
      mutate(group_value = clean_plot_label(.data[[group_col]])) %>%
      count(pfam_id, group_value, name = "n") %>%
      group_by(pfam_id) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
      mutate(pfam_id = factor(pfam_id, levels = pfam_levels))
    p <- ggplot(d, aes(x = pfam_id, y = prop, fill = group_value)) +
      geom_col() + coord_flip() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = NULL, y = "Share within PFAM", fill = clean_plot_label(group_col), title = title) +
      theme_minimal(base_size = 10)
    save_plot_pair(p, stem, output_dir, 9.5, 7.2, formats)
  }
  stacked_plot("phylum", "top_pfam_domains_overlapping_lcr_by_phylum", "PFAM-LCR overlap domains split by phylum")
  stacked_plot("lifestyle", "top_pfam_domains_overlapping_lcr_by_lifestyle", "PFAM-LCR overlap domains split by lifestyle")
}

if ("lcr_total_length" %in% names(plot_data)) {
  df_lcr <- plot_data %>% filter(is.finite(lcr_total_length), is.finite(tAI_z), !is.na(overlap_label))
  if (nrow(df_lcr) >= 50L) {
    ann <- lm_annotation(df_lcr$lcr_total_length, df_lcr$tAI_z)
    p_lcr <- sample_rows(df_lcr, max_plot_rows) %>%
      ggplot(aes(x = lcr_total_length, y = tAI_z, colour = overlap_label)) +
      geom_point(alpha = 0.12, size = 0.45) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
      annotate_top_right(ann) +
      labs(
        x = "Total LCR length [aa]", y = "within-genome standardized tAI (tAI_z)", colour = NULL,
        title = "LCR length vs tAI by PFAM-LCR overlap status",
        subtitle = "TRUE/FALSE labels replaced by explicit biological labels. Regression annotation is global."
      ) + theme_minimal(base_size = 11)
    save_plot_pair(p_lcr, "lcr_length_vs_tai_z_by_pfam_overlap", output_dir, 8.2, 5.8, formats)
  }
}

readr::write_tsv(
  plot_data %>% count(structural_category, name = "n_genes"),
  file.path(output_dir, "pfam_lcr_structural_category_summary.tsv")
)
message("Done: ", output_dir)
