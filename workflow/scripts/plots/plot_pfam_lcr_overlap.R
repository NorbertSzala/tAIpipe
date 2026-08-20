#!/usr/bin/env Rscript
# PFAM-LCR overlap plots.
# tAI_z is the default for cross-genome plots because it removes genome-specific
# shifts in tAI scale. Plain tAI plots are also created when tAI is available.

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
    pfam_tail_enrichment = "results/statistics/pfam_tai_tail_enrichment.tsv",
    top_n = "25",
    min_group_n = "5",
    formats = "png,pdf",
    max_plot_rows = "160000",
    seed = "1",
    manifest_output = "results/plots/pfam_lcr_overlap/plot_manifest.tsv"
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
    pfam_tail_enrichment = snakemake@input[["pfam_tail_enrichment"]] %||% "results/statistics/pfam_tai_tail_enrichment.tsv",
    output_dir = snakemake@params[["output_dir"]] %||% "results/plots/pfam_lcr_overlap",
    pfam_description_table = snakemake@params[["pfam_description_table"]] %||% "",
    top_n = as.character(snakemake@params[["top_n"]] %||% snakemake@params[["top_n_pfam"]] %||% 25),
    min_group_n = as.character(snakemake@params[["min_group_n"]] %||% 5),
    formats = paste(snakemake@params[["formats"]] %||% c("png", "pdf"), collapse = ","),
    max_plot_rows = as.character(snakemake@params[["max_plot_rows"]] %||% 160000),
    seed = as.character(snakemake@params[["seed"]] %||% 1)
  )
} else parse_cli()

output_dir <- args$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(output_dir, c(
  "pfam_lcr_structural_categories_tai_z_boxplot.png",
  "pfam_lcr_structural_categories_tai_z_boxplot.pdf",
  "pfam_lcr_overlap_tai_z_violin.png", "pfam_lcr_overlap_tai_z_violin.pdf",
  "pfam_lcr_overlap_raw_tai_violin.png", "pfam_lcr_overlap_raw_tai_violin.pdf",
  "binary_feature_split_violins_tai_z.png", "binary_feature_split_violins_tai_z.pdf",
  "binary_feature_split_violins_tai.png", "binary_feature_split_violins_tai.pdf"
)), force = TRUE)
formats <- strsplit(args$formats, ",", fixed = TRUE)[[1]] |> trimws()
top_n <- as.integer(args$top_n_pfam %||% args$top_n)
min_group_n <- as.integer(args$min_group_n)
max_plot_rows <- as.integer(args$max_plot_rows)
set.seed(as.integer(args$seed))

sample_rows <- function(df, max_n) {
  if (!is.finite(max_n) || max_n <= 0L) return(df)
  k <- min(nrow(df), max_n)
  if (k <= 0L) return(df[0, , drop = FALSE])
  if (nrow(df) <= k) return(df)
  dplyr::slice_sample(df, n = k, replace = FALSE)
}

write_plot_manifest <- function() {
  manifest_files <- list.files(
    output_dir,
    pattern = "\\.(png|pdf|tsv|txt)$",
    full.names = TRUE
  )
  manifest <- tibble::tibble(
    file = basename(manifest_files),
    size_bytes = file.info(manifest_files)$size
  ) %>% arrange(file)
  readr::write_tsv(manifest, args$manifest_output)
}

write_missing <- function(text) {
  writeLines(
    text,
    file.path(output_dir, "XXXXX_PFAM_LCR_PLOTS_NOT_CREATED.txt")
  )
  write_plot_manifest()
  message(text)
}

term_present <- function(x) {
  y <- trimws(as.character(x))
  !is.na(y) & nzchar(y) & !tolower(y) %in% c("na", "none", "null", "false", "0")
}

normalize_pfam_id <- function(x) {
  stringr::str_extract(toupper(as.character(x)), "PF[0-9]{5}")
}

read_pfam_descriptions <- function(path) {
  out <- tibble(pfam_id = character(), pfam_description = character())
  if (!nzchar(path) || !file.exists(path)) return(out)

  # Supported inputs:
  #   1) delimited table (pfam_id/accession + description/name/function),
  #   2) Pfam-A.hmm records (NAME/ACC/DESC),
  #   3) Pfam-A.hmm.dat Stockholm records (#=GF ID/AC/DE).
  first_lines <- readLines(path, n = 80L, warn = FALSE)
  is_hmm <- any(grepl("^ACC\\s+PF[0-9]{5}", first_lines))
  is_hmm_dat <- any(grepl("^(#=GF\\s+)?AC\\s+PF[0-9]{5}", first_lines))

  if (is_hmm || is_hmm_dat) {
    lines <- readLines(path, warn = FALSE)
    records <- strsplit(paste(lines, collapse = "\n"), "\n//", fixed = FALSE)[[1]]
    parsed <- lapply(records, function(rec) {
      ls <- strsplit(rec, "\n", fixed = TRUE)[[1]]
      accession_line <- grep("^(#=GF\\s+)?(ACC|AC)\\s+PF[0-9]{5}", ls, value = TRUE)
      description_line <- grep("^(#=GF\\s+)?(DESC|DE)\\s+", ls, value = TRUE)
      name_line <- grep("^(#=GF\\s+)?(NAME|ID)\\s+", ls, value = TRUE)
      accession <- if (length(accession_line) > 0L) sub("^(#=GF\\s+)?(ACC|AC)\\s+", "", accession_line[[1]]) else NA_character_
      description <- if (length(description_line) > 0L) sub("^(#=GF\\s+)?(DESC|DE)\\s+", "", description_line[[1]]) else NA_character_
      name <- if (length(name_line) > 0L) sub("^(#=GF\\s+)?(NAME|ID)\\s+", "", name_line[[1]]) else NA_character_
      tibble(
        pfam_id = normalize_pfam_id(accession),
        pfam_description = stringr::str_squish(dplyr::coalesce(description, name))
      )
    })
    return(
      bind_rows(parsed) %>%
        filter(!is.na(pfam_id), !is.na(pfam_description), nzchar(pfam_description)) %>%
        distinct(pfam_id, .keep_all = TRUE)
    )
  }

  tab <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  id_col <- intersect(c("pfam_id", "accession", "pfam_acc", "AC", "ACC"), names(tab))
  if (length(id_col) == 0L) stop("PFAM description table must contain pfam_id/accession")
  desc_col <- intersect(c("pfam_description", "description", "pfam_name", "name", "function", "DE", "DESC", "NAME"), names(tab))
  if (length(desc_col) == 0L) stop("PFAM description table must contain a description/name/function column")
  tab %>%
    transmute(
      pfam_id = normalize_pfam_id(.data[[id_col[[1]]]]),
      pfam_description = stringr::str_squish(stringr::str_replace_all(as.character(.data[[desc_col[[1]]]]), "[\r\n\t]+", " "))
    ) %>%
    filter(!is.na(pfam_id), nzchar(pfam_description), !tolower(pfam_description) %in% c("na", "none", "description not supplied")) %>%
    distinct(pfam_id, .keep_all = TRUE)
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

pfam_present_coerced <- coerce_binary_plot(genes$pfam_present)
lcr_present_coerced <- coerce_binary_plot(genes$lcr_present)

plot_data <- genes %>%
  mutate(
    tAI_z = suppressWarnings(as.numeric(tAI_z)),
    tAI = suppressWarnings(as.numeric(tAI)),
    pfam_lcr_overlap_present = coerce_binary_plot(pfam_lcr_overlap_present),
    pfam_present_bool = ifelse(!is.na(pfam_present_coerced), pfam_present_coerced, term_present(pfam_terms)),
    lcr_present_bool = ifelse(!is.na(lcr_present_coerced), lcr_present_coerced, is.finite(suppressWarnings(as.numeric(lcr_total_length))) & suppressWarnings(as.numeric(lcr_total_length)) > 0),
    lcr_total_length = suppressWarnings(as.numeric(lcr_total_length)),
    phylum = clean_plot_label(phylum),
    lifestyle = clean_plot_label(lifestyle),
    overlap_label = factor(
      if_else(pfam_lcr_overlap_present, "PFAM-LCR overlap", "no PFAM-LCR overlap"),
      levels = c("no PFAM-LCR overlap", "PFAM-LCR overlap")
    ),
    structural_category = case_when(
      pfam_lcr_overlap_present ~ "PFAM-LCR overlap",
      pfam_present_bool & lcr_present_bool ~ "PFAM and LCR, no overlap",
      pfam_present_bool & !lcr_present_bool ~ "PFAM only",
      !pfam_present_bool & lcr_present_bool ~ "LCR only",
      TRUE ~ "Neither PFAM nor LCR"
    )
  )

if (sum(plot_data$pfam_lcr_overlap_present, na.rm = TRUE) < 10L) {
  write_missing("XXXXX: fewer than 10 genes have PFAM-LCR overlap; plot would not be stable.")
  quit(save = "no", status = 0)
}

plot_binary_feature_grid <- function() {
  feature_map <- c(
    "PFAM domains" = "pfam_present_bool",
    "LCR" = "lcr_present_bool",
    "PFAM-LCR overlap" = "pfam_lcr_overlap_present"
  )
  metric_map <- c("tAI z-score" = "tAI_z", "Raw tAI" = "tAI")
  metric_map <- metric_map[metric_map %in% names(plot_data)]
  d <- purrr::imap_dfr(metric_map, function(metric_col, metric_label) {
    purrr::imap_dfr(feature_map, function(feature_col, feature_label) {
      plot_data %>%
        transmute(
          sample,
          metric = metric_label,
          feature = feature_label,
          status = .data[[feature_col]],
          value = suppressWarnings(as.numeric(.data[[metric_col]]))
        )
    })
  }) %>%
    filter(!is.na(status), is.finite(value)) %>%
    mutate(
      status = factor(if_else(status, "Present", "Absent"), levels = c("Absent", "Present")),
      split_side = split_side_from_level(status, "Absent"),
      metric = factor(metric, levels = names(metric_map)),
      feature = factor(feature, levels = names(feature_map))
    )
  if (nrow(d) < 20L) return(invisible(NULL))
  counts <- d %>%
    count(metric, feature, status, name = "n") %>%
    mutate(
      count_x = if_else(status == "Absent", 0.82, 1.18),
      count_label = format(n, scientific = FALSE, trim = TRUE, big.mark = "")
    )
  d <- sample_rows(d, max_plot_rows)
  p <- ggplot(
    d,
    aes(
      x = 1, y = value, fill = status,
      group = interaction(metric, feature, status), split_side = split_side
    )
  ) +
    geom_split_violin_project(trim = TRUE, alpha = 0.72, width = 0.92, scale = "width", colour = "grey30", linewidth = 0.22) +
    geom_boxplot(
      aes(group = interaction(metric, feature, status)),
      width = 0.075,
      position = position_dodge(width = 0.18),
      outlier.shape = NA,
      fill = "white",
      colour = "grey25",
      alpha = 0.92,
      linewidth = 0.35
    ) +
    geom_text(
      data = counts,
      aes(x = count_x, y = -Inf, label = count_label),
      inherit.aes = FALSE,
      vjust = 2.1,
      size = 3.5,
      colour = "grey25"
    ) +
    facet_grid(metric ~ feature, scales = "free_y") +
    scale_fill_manual(values = binary_grey_values(c("Absent", "Present")), drop = FALSE) +
    scale_x_continuous(breaks = NULL) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = "Metric value",
      fill = NULL,
      title = "tAI distributions by binary structural annotations",
      subtitle = "Numbers below each violin half are gene counts before plotting downsampling."
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      strip.text = element_text(face = "bold", size = 12.5), legend.position = "bottom",
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12.5),
      panel.spacing = grid::unit(1.1, "lines"), plot.margin = margin(8, 8, 28, 8)
    )
  save_plot_pair(p, "binary_feature_split_violins_tai_and_tai_z", output_dir, 14.8, 10.2, formats)
}
plot_binary_feature_grid()

# Structural categories: use violin, not boxplot only. Wrap long labels.
df_cat <- plot_data %>% filter(is.finite(tAI_z), !is.na(structural_category))
p_cat <- safe_kruskal_p(df_cat$tAI_z, df_cat$structural_category)
cat_plot <- df_cat %>%
  mutate(structural_category = stringr::str_wrap(structural_category, width = 22)) %>%
  sample_rows(max_plot_rows) %>%
  ggplot(aes(x = structural_category, y = tAI_z, fill = structural_category)) +
  geom_violin(trim = TRUE, alpha = 0.68, width = 1.12, scale = "width") +
  geom_boxplot(width = 0.13, outlier.alpha = 0.18, alpha = 0.90) +
  annotate_top_right(sig_label(p_cat)) +
  labs(
    x = NULL,
    y = "tAI z-score",
    title = "tAI z-score across PFAM/LCR structural annotation categories"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11.5, lineheight = 0.95),
    axis.text.y = element_text(size = 11.5),
    plot.title = element_text(size = 16, face = "bold")
  )
save_plot_pair(cat_plot, "pfam_lcr_structural_categories_tai_z_violin", output_dir, 10.5, 6.2, formats)

# Long PFAM table.
pfam_long <- plot_data %>%
  filter(pfam_lcr_overlap_present, term_present(pfam_lcr_overlap_terms)) %>%
  transmute(sample, gene_id, phylum, lifestyle, tAI_z, pfam_token = str_split(as.character(pfam_lcr_overlap_terms), ";|\\|")) %>%
  tidyr::unnest_longer(pfam_token) %>%
  mutate(
    pfam_token = str_squish(pfam_token),
    pfam_id = normalize_pfam_id(pfam_token),
    token_description = str_squish(str_remove(pfam_token, regex("PF[0-9]{5}(\\.[0-9]+)?", ignore_case = TRUE))),
    token_description = str_remove(token_description, "^[[:punct:]\\s]+"),
    token_description = na_if(token_description, "")
  ) %>%
  filter(!is.na(pfam_id)) %>%
  distinct(sample, gene_id, pfam_id, .keep_all = TRUE)

pfam_desc <- read_pfam_descriptions(args$pfam_description_table)

if (nrow(pfam_long) > 0L) {
  pfam_summary <- pfam_long %>%
    group_by(pfam_id) %>%
    summarise(
      n_genes = n_distinct(paste(sample, gene_id, sep = "\t")),
      n_samples = n_distinct(sample),
      median_tAI_z = median(tAI_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      pfam_long %>% group_by(pfam_id) %>% summarise(token_description = first(na.omit(token_description), default = NA_character_), .groups = "drop"),
      by = "pfam_id"
    ) %>%
    left_join(pfam_desc, by = "pfam_id") %>%
    mutate(pfam_description = dplyr::coalesce(pfam_description, token_description, "Description not supplied")) %>%
    arrange(desc(n_genes), desc(n_samples))

  readr::write_tsv(pfam_summary, file.path(output_dir, "pfam_lcr_overlap_domain_summary.tsv"), na = "NA")
  readr::write_tsv(
    pfam_summary %>% filter(pfam_description == "Description not supplied") %>% select(pfam_id, n_genes, n_samples),
    file.path(output_dir, "pfam_lcr_overlap_missing_pfam_descriptions.tsv"),
    na = "NA"
  )

  top_ids <- head(pfam_summary$pfam_id, top_n)
  base_counts <- pfam_long %>% filter(pfam_id %in% top_ids)
  pfam_levels <- pfam_summary %>% filter(pfam_id %in% top_ids) %>% arrange(n_genes) %>% pull(pfam_id)

  label_map <- pfam_summary %>%
    mutate(label = paste0(pfam_id, " - ", pfam_description)) %>%
    select(pfam_id, label)

  total_plot <- pfam_summary %>%
    filter(pfam_id %in% top_ids) %>%
    mutate(pfam_id = factor(pfam_id, levels = pfam_levels)) %>%
    ggplot(aes(x = pfam_id, y = n_genes)) +
    geom_col(fill = "grey55", colour = "grey30", linewidth = 0.15, width = 0.86) +
    coord_flip() +
    scale_x_discrete(labels = setNames(stringr::str_wrap(label_map$label, 58), label_map$pfam_id)) +
    labs(x = NULL, y = "Genes with PFAM-LCR overlap", title = "Most frequent PFAM domains overlapping LCRs") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 12.5, lineheight = 1.05),
      axis.text.x = element_text(size = 12),
      plot.title = element_text(size = 17, face = "bold")
    )
  save_plot_pair(
    total_plot,
    "top_pfam_domains_overlapping_lcr",
    output_dir,
    12.5,
    max(10.0, 3.8 + 0.58 * length(top_ids)),
    formats
  )

  prevalence_plot <- function(group_col, stem, title) {
    # Normalize for unequal clade/lifestyle sizes: each value is the percentage
    # of genomes in a group containing >=1 gene with the PFAM-LCR overlap.
    all_meta <- plot_data %>%
      distinct(sample, group_value = .data[[group_col]]) %>%
      mutate(group_value = clean_plot_label(group_value))
    raw_sizes <- all_meta %>% count(group_value, name = "n_genomes")
    keep_groups <- raw_sizes %>% filter(n_genomes >= min_group_n) %>% pull(group_value)
    all_meta <- all_meta %>%
      mutate(group_value = if_else(group_value %in% keep_groups, group_value, "Other"))

    denominators <- all_meta %>% count(group_value, name = "n_genomes")
    present <- base_counts %>%
      transmute(sample, pfam_id, group_value = clean_plot_label(.data[[group_col]])) %>%
      mutate(group_value = if_else(group_value %in% keep_groups, group_value, "Other")) %>%
      distinct(sample, pfam_id, group_value) %>%
      count(pfam_id, group_value, name = "n_genomes_with_domain")

    d <- tidyr::crossing(pfam_id = top_ids, group_value = denominators$group_value) %>%
      left_join(present, by = c("pfam_id", "group_value")) %>%
      left_join(denominators, by = "group_value") %>%
      mutate(
        n_genomes_with_domain = coalesce(n_genomes_with_domain, 0L),
        prevalence_pct = 100 * n_genomes_with_domain / n_genomes
      )

    group_order <- order_project_groups(denominators$group_value, group_col)
    group_info <- denominators %>%
      mutate(
        group_value = factor(group_value, levels = group_order),
        group_label = paste0(as.character(group_value), " (", n_genomes, ")")
      ) %>% arrange(group_value)
    label_lookup <- stats::setNames(group_info$group_label, as.character(group_info$group_value))
    d <- d %>%
      mutate(
        pfam_id = factor(pfam_id, levels = pfam_levels),
        group_value = factor(group_value, levels = group_order),
        group_label = factor(label_lookup[as.character(group_value)], levels = group_info$group_label)
      )

    colour_values <- project_category_colours(as.character(group_info$group_value))
    names(colour_values) <- group_info$group_label
    p <- ggplot(d, aes(x = pfam_id, y = prevalence_pct, fill = group_label)) +
      geom_col(position = position_dodge2(width = 0.98, preserve = "single"), width = 0.88,
               colour = "grey25", linewidth = 0.10) +
      coord_flip() +
      scale_x_discrete(labels = setNames(stringr::str_wrap(label_map$label, 52), label_map$pfam_id)) +
      scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100), expand = expansion(mult = c(0, 0.03))) +
      scale_fill_manual(values = colour_values, drop = FALSE) +
      labs(
        x = NULL, y = "Genomes containing the overlap domain [%]",
        fill = paste0(clean_plot_label(group_col), " (number of genomes)"),
        title = title,
        subtitle = "Prevalence is normalized within each group; one genome contributes at most once per PFAM domain."
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "vertical", legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 11.5),
        axis.text.y = element_text(size = 12.5, lineheight = 1.05),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 17, face = "bold"),
        plot.subtitle = element_text(size = 11.5)
      ) +
      guides(fill = guide_legend(nrow = 3, byrow = TRUE))
    save_plot_pair(p, stem, output_dir, 15.3, max(12.0, 4.2 + 0.62 * length(top_ids)), formats)

    # Complementary view: retain the absolute number of genomes as bar length,
    # but encode within-group prevalence by fill. Faceting keeps group identity
    # explicit and avoids comparing unlabeled dodge positions.
    p_counts <- ggplot(
      d,
      aes(x = pfam_id, y = n_genomes_with_domain, fill = prevalence_pct)
    ) +
      geom_col(width = 0.96, colour = "grey25", linewidth = 0.12) +
      coord_flip() +
      facet_wrap(vars(group_label), scales = "free_y", ncol = 2) +
      scale_x_discrete(
        limits = pfam_levels,
        drop = FALSE,
        labels = setNames(
          stringr::str_wrap(label_map$label, 42),
          label_map$pfam_id
        )
      ) +
      scale_fill_viridis_c(
        option = "C",
        begin = 0.12,
        end = 0.88,
        limits = c(0, 100),
        name = "Within-group prevalence [%]",
        guide = guide_colourbar(
          title.position = "top",
          barwidth = grid::unit(36, "cm"),
          barheight = grid::unit(0.75, "cm")
        )
      ) +
      labs(
        x = NULL,
        y = "Genomes containing the overlap domain",
        title = stringr::str_wrap(
          paste0(title, ": counts coloured by normalized prevalence"),
          width = 64
        )
      ) +
      theme_minimal(base_size = 15) +
      theme(
        axis.text.y = element_text(size = 14.5, lineheight = 1.08),
        axis.text.x = element_text(size = 13.0),
        strip.text = element_text(face = "bold", size = 13.2),
        legend.position = "bottom",
        legend.title = element_text(size = 13.5, face = "bold"),
        legend.text = element_text(size = 13.0),
        plot.title = element_text(size = 18, face = "bold", lineheight = 1.05)
      )
    save_plot_pair(
      p_counts,
      paste0(stem, "_counts_coloured_by_prevalence"),
      output_dir,
      18.5,
      max(14.0, 4.5 + 0.62 * length(top_ids)),
      formats
    )

    readr::write_tsv(
      d %>% transmute(pfam_id = as.character(pfam_id), group = as.character(group_value),
                      n_genomes_with_domain, n_genomes, prevalence_pct),
      file.path(output_dir, paste0(stem, "_normalized_prevalence.tsv")), na = "NA"
    )
  }

  prevalence_plot("phylum", "top_pfam_domains_overlapping_lcr_by_phylum", "PFAM-LCR overlap-domain prevalence by phylum")
  prevalence_plot("lifestyle", "top_pfam_domains_overlapping_lcr_by_lifestyle", "PFAM-LCR overlap-domain prevalence by lifestyle")
}

# PFAM enrichment in high/low tAI tails is computed separately with genomes as
# CMH strata. This panel answers a different question than PFAM-LCR overlap: it
# compares each within-genome tail with the remaining genes in that genome.
if (file.exists(args$pfam_tail_enrichment)) {
  pfam_tail <- readr::read_tsv(args$pfam_tail_enrichment, show_col_types = FALSE, progress = FALSE)
  required_tail <- c(
    "PFAM", "description", "n_count", "mean_tAI", "label",
    "tail_percent", "log2_enrichment", "q_value"
  )
  if (all(required_tail %in% names(pfam_tail)) && nrow(pfam_tail) > 0L) {
    pfam_tail <- pfam_tail %>%
      mutate(
        n_count = suppressWarnings(as.numeric(n_count)),
        mean_tAI = suppressWarnings(as.numeric(mean_tAI)),
        tail_percent = suppressWarnings(as.numeric(tail_percent)),
        log2_enrichment = suppressWarnings(as.numeric(log2_enrichment)),
        q_value = suppressWarnings(as.numeric(q_value)),
        label = factor(label, levels = c("highest", "lowest")),
        tail_label = paste0(format(tail_percent, trim = TRUE), "% tail"),
        tail_label = factor(tail_label, levels = c("10% tail", "1% tail")),
        domain_label = paste0(PFAM, " - ", description)
      ) %>%
      filter(is.finite(n_count), n_count > 0, is.finite(mean_tAI), !is.na(label), !is.na(tail_label)) %>%
      group_by(label, tail_label) %>%
      arrange(q_value, desc(abs(log2_enrichment)), desc(n_count), .by_group = TRUE) %>%
      slice_head(n = top_n) %>%
      ungroup()

    finite_effect <- abs(pfam_tail$log2_enrichment[is.finite(pfam_tail$log2_enrichment)])
    display_limit <- if (length(finite_effect) > 0L) {
      min(8, max(2, as.numeric(stats::quantile(finite_effect, 0.98, na.rm = TRUE, names = FALSE))))
    } else {
      2
    }
    pfam_tail <- pfam_tail %>%
      mutate(
        plot_log2_enrichment = pmax(-display_limit, pmin(display_limit, log2_enrichment)),
        clipped = !is.finite(log2_enrichment) | abs(log2_enrichment) > display_limit,
        domain_label = factor(
          domain_label,
          levels = unique(domain_label[order(plot_log2_enrichment, na.last = TRUE)])
        )
      )

    readr::write_tsv(
      pfam_tail %>% mutate(domain_label = as.character(domain_label)),
      file.path(output_dir, "pfam_tai_tail_enrichment_domains_used_for_plot.tsv"),
      na = "NA"
    )

    p_pfam_tail <- ggplot(
      pfam_tail,
      aes(
        x = plot_log2_enrichment, y = domain_label,
        size = n_count, fill = mean_tAI
      )
    ) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.45) +
      geom_point(shape = 21, alpha = 0.94, colour = "grey15", stroke = 0.45) +
      facet_grid(label ~ tail_label, scales = "free_y", space = "free_y") +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, 52)) +
      scale_fill_viridis_c(
        option = "inferno", begin = 0.08, end = 0.92,
        name = "Mean raw tAI",
        guide = guide_colourbar(
          title.position = "top",
          barwidth = grid::unit(13, "cm"),
          barheight = grid::unit(0.70, "cm")
        )
      ) +
      scale_size_continuous(range = c(3.0, 10.5), name = "Tail genes with PFAM") +
      labs(
        x = "log2(CMH common odds ratio)", y = NULL,
        title = "PFAM domains enriched or depleted in within-genome tAI tails",
        subtitle = paste0(
          "Top and bottom tails are selected separately in every genome. Effects beyond +/-",
          formatC(display_limit, digits = 2, format = "fg"), " are clipped for display only."
        )
      ) +
      theme_minimal(base_size = 15) +
      theme(
        axis.text.y = element_text(size = 13.3, lineheight = 1.03),
        axis.text.x = element_text(size = 13.5),
        axis.title = element_text(size = 14.5),
        strip.text = element_text(size = 14.0, face = "bold"),
        strip.background = element_rect(fill = "grey94", colour = "grey60", linewidth = 0.3),
        legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size = 14.5, face = "bold"),
        legend.text = element_text(size = 14.0),
        plot.title = element_text(size = 19, face = "bold"),
        plot.subtitle = element_text(size = 12.8),
        panel.spacing = grid::unit(1.2, "lines")
      )
    save_plot_pair(
      p_pfam_tail,
      "pfam_tai_tail_enrichment_dotplot",
      output_dir,
      18.0,
      max(12.0, 5.0 + 0.34 * nrow(pfam_tail)),
      formats
    )
  }
}

if ("lcr_total_length" %in% names(plot_data)) {
  df_lcr_all <- plot_data %>%
    filter(
      is.finite(lcr_total_length),
      lcr_total_length > 0,
      is.finite(tAI_z),
      !is.na(overlap_label)
    )

  plot_lcr_scatter <- function(df_lcr, stem, title_suffix = NULL) {
    if (nrow(df_lcr) < 50L) return(invisible(NULL))

    ann_tbl <- df_lcr %>%
      group_by(overlap_label) %>%
      summarise(label = lm_annotation(lcr_total_length, tAI_z), .groups = "drop") %>%
      mutate(label = paste0(as.character(overlap_label), ": ", label))
    ann <- paste(ann_tbl$label, collapse = "\n")

    point_colours <- binary_grey_values(levels(df_lcr$overlap_label))
    scatter_counts <- df_lcr %>% count(overlap_label, name = "n")
    scatter_labels <- stats::setNames(
      paste0(scatter_counts$overlap_label, " (", format(scatter_counts$n, big.mark = ","), ")"),
      scatter_counts$overlap_label
    )

    p_lcr <- sample_rows(df_lcr, max_plot_rows) %>%
      ggplot(aes(x = lcr_total_length, y = tAI_z, colour = overlap_label)) +
      geom_point(alpha = 0.34, size = 0.55) +
      geom_smooth(
        aes(group = overlap_label, colour = overlap_label),
        method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.75
      ) +
      annotate_top_right(ann, size = 3.05) +
      scale_colour_manual(values = point_colours, labels = scatter_labels, drop = FALSE) +
      labs(
        x = "Total LCR length [aa]",
        y = "tAI z-score",
        colour = NULL,
        title = if (is.null(title_suffix) || !nzchar(title_suffix)) {
          "LCR length vs tAI z-score by PFAM-LCR overlap status"
        } else {
          "LCR length vs tAI z-score by PFAM-LCR overlap status; lowest 90% of lengths"
        }
      ) +
      theme_minimal(base_size = 15) +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14.2),
        legend.text = element_text(size = 13),
        legend.key.width = grid::unit(0.9, "cm"),
        legend.spacing.x = grid::unit(0.40, "cm"),
        plot.title = element_text(size = 18, face = "bold")
      )

    p_lcr_full <- add_marginal_densities(
      p_lcr, df_lcr, lcr_total_length, tAI_z, overlap_label, point_colours,
      top_height = 1.4, right_width = 1.5
    )
    save_plot_pair(p_lcr_full, stem, output_dir, 12.2, 8.4, formats)
  }

  plot_lcr_scatter(
    df_lcr_all,
    "lcr_length_vs_tai_z_by_pfam_overlap"
  )

  if (nrow(df_lcr_all) >= 50L) {
    q90 <- as.numeric(stats::quantile(df_lcr_all$lcr_total_length, 0.90, na.rm = TRUE, names = FALSE))
    df_lcr_90 <- df_lcr_all %>% filter(lcr_total_length <= q90)
    readr::write_tsv(
      tibble::tibble(
        filter = "lowest 90% of positive total LCR lengths",
        quantile = 0.90,
        cutoff_aa = q90,
        rows_before = nrow(df_lcr_all),
        rows_after = nrow(df_lcr_90)
      ),
      file.path(output_dir, "lcr_length_vs_tai_z_lowest90_filter.tsv")
    )
    plot_lcr_scatter(
      df_lcr_90,
      "lcr_length_vs_tai_z_by_pfam_overlap_lowest90pct",
      " (lowest 90% of LCR lengths)"
    )
  }
}

readr::write_tsv(
  plot_data %>% count(structural_category, name = "n_genes"),
  file.path(output_dir, "pfam_lcr_structural_category_summary.tsv"),
  na = "NA"
)

write_plot_manifest()
message("Done: ", output_dir)
