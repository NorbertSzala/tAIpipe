#!/usr/bin/env Rscript
# Genome-level GC/GC3s/tAI plots. Uses genomes as independent plotting units.
# This version avoids ggplot2 position_dodge warnings by not mapping fill to the
# internal boxplot layer. Violin widths are wider and subtitles are wrapped.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})
source("workflow/scripts/lib/plot_style_helpers.R")

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    genome_summary = "results/tables/genome_summary.tsv",
    output_dir = "results/plots/script_suggestions/genome_gc_tai",
    min_n = "5",
    formats = "png,pdf"
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
    genome_summary = snakemake@input[["genome_summary"]],
    output_dir = snakemake@params[["output_dir"]] %||% "results/plots/script_suggestions/genome_gc_tai",
    min_n = as.character(snakemake@params[["min_n"]] %||% 5),
    formats = paste(snakemake@params[["formats"]] %||% c("png", "pdf"), collapse = ",")
  )
} else parse_args()

formats <- strsplit(args$formats, ",", fixed = TRUE)[[1]] |> trimws()
min_n <- as.integer(args$min_n)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
missing <- tibble(plot = character(), reason = character())
add_missing <- function(plot, reason) {
  missing <<- bind_rows(missing, tibble(plot = plot, reason = reason))
  message("XXXXX skipped ", plot, ": ", reason)
}

if (!file.exists(args$genome_summary)) stop("Missing genome summary: ", args$genome_summary)
genomes <- readr::read_tsv(args$genome_summary, show_col_types = FALSE, progress = FALSE)
if (!"phylum" %in% names(genomes)) genomes$phylum <- NA_character_
if (!"lifestyle" %in% names(genomes)) genomes$lifestyle <- NA_character_

genomes <- genomes %>%
  mutate(
    phylum = as.character(phylum),
    lifestyle = as.character(lifestyle),
    dikarya_status = case_when(
      phylum %in% c("Ascomycota", "Basidiomycota") ~ "Dikarya",
      is.na(phylum) | !nzchar(phylum) ~ NA_character_,
      TRUE ~ "non-Dikarya"
    )
  )

plot_violin_group <- function(metric, group, stem, y_label, collapse = FALSE) {
  if (!all(c(metric, group) %in% names(genomes))) {
    add_missing(stem, paste("missing", metric, "or", group)); return(NULL)
  }

  df <- genomes %>%
    transmute(group = .data[[group]], value = suppressWarnings(as.numeric(.data[[metric]]))) %>%
    filter(!is.na(group), nzchar(group), is.finite(value))

  if (collapse) df <- collapse_small_groups(df, "group", min_n = min_n, other_label = "Other")
  if (nrow(df) < 3L || n_distinct(df$group) < 2L) {
    add_missing(stem, "too few groups"); return(NULL)
  }

  other <- attr(df, "other_members")
  labs_map <- add_group_n(df, "group")$labels
  pval <- if (n_distinct(df$group) == 2L) safe_wilcox_p(df$value, df$group) else safe_kruskal_p(df$value, df$group)
  lev <- df %>% group_by(group) %>% summarise(m = median(value), .groups = "drop") %>% arrange(m) %>% pull(group)
  df$group <- factor(df$group, levels = lev)

  subtitle <- if (!is.null(other) && length(other) > 0L) {
    wrap_text(paste0("Other contains: ", paste(clean_plot_label(other), collapse = ", ")), 85)
  } else {
    wrap_text("Each point is one genome; categories are ordered by median.", 85)
  }

  cols <- project_category_colours(levels(df$group))

  p <- ggplot(df, aes(x = group, y = value)) +
    geom_violin(aes(fill = group), trim = FALSE, alpha = 0.70, width = 1.12, position = position_identity()) +
    geom_boxplot(width = 0.07, outlier.alpha = 0.35, alpha = 0.90,
                 fill = "white", colour = "#222222", position = position_identity()) +
    annotate_top_right(sig_label(pval)) +
    scale_x_discrete(labels = labs_map) +
    scale_fill_manual(values = cols, guide = "none") +
    labs(x = NULL, y = y_label, title = paste(y_label, "by", clean_plot_label(group)), subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.subtitle = element_text(lineheight = 0.95))

  save_plot_pair(p, stem, args$output_dir, 10, 6, formats)
  p + labs(title = y_label)
}

p1 <- plot_violin_group("genome_gc", "phylum", "genome_GC_by_phylum", "Genome GC fraction", collapse = TRUE)
p2 <- plot_violin_group("mean_GC3s", "phylum", "CDS_GC3s_by_phylum", "Mean CDS GC3s", collapse = TRUE)
p3 <- plot_violin_group("mean_tAI", "dikarya_status", "dikarya_non_dikarya_tAI_CDS", "Mean tAI", collapse = FALSE)
p4 <- plot_violin_group("mean_GC3s", "dikarya_status", "dikarya_non_dikarya_GC3s_CDS", "Mean CDS GC3s", collapse = FALSE)
p5 <- plot_violin_group("genome_gc", "dikarya_status", "dikarya_non_dikarya_genome_GC", "Genome GC fraction", collapse = FALSE)

if (requireNamespace("patchwork", quietly = TRUE)) {
  panels <- list(p1, p2, p3, p4, p5)
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (length(panels) > 1L) {
    grid <- patchwork::wrap_plots(panels, ncol = 2)
    save_plot_pair(grid, "genome_metric_violins_grid", args$output_dir, 13, 9, formats)
  }
}

if (all(c("genome_gc", "mean_GC3s") %in% names(genomes))) {
  df <- genomes %>%
    transmute(
      genome_gc = as.numeric(genome_gc),
      mean_GC3s = as.numeric(mean_GC3s),
      dikarya_status = dikarya_status
    ) %>%
    filter(is.finite(genome_gc), is.finite(mean_GC3s))

  if (nrow(df) >= 3L) {
    ann <- lm_annotation(df$genome_gc, df$mean_GC3s)
    p <- ggplot(df, aes(x = genome_gc, y = mean_GC3s, colour = dikarya_status)) +
      geom_point(size = 2, alpha = 0.8) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
      annotate_top_right(ann) +
      labs(x = "Genome GC fraction", y = "Mean CDS GC3s", colour = NULL, title = "Genome GC vs CDS GC3s") +
      theme_minimal(base_size = 11)
    save_plot_pair(p, "genome_GC_vs_CDS_GC3s", args$output_dir, 8, 5.7, formats)
  }
}

if (nrow(missing) == 0L) missing <- tibble(plot = "none", reason = "all possible genome GC/tAI plots generated")
readr::write_tsv(missing, file.path(args$output_dir, "XXXXX_missing_plots.tsv"))
message("Done: ", args$output_dir)
