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
unlink(file.path(args$output_dir, c(
  "dikarya_non_dikarya_GC3s_CDS.png", "dikarya_non_dikarya_GC3s_CDS.pdf",
  "dikarya_non_dikarya_genome_GC.png", "dikarya_non_dikarya_genome_GC.pdf"
)), force = TRUE)
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
    wrap_text("Each observation is one genome; categories are ordered by median.", 85)
  }

  cols <- if (identical(group, "dikarya_status")) {
    c("non-Dikarya" = "#B8B8B8", "Dikarya" = "#1F4E79")[levels(df$group)]
  } else {
    project_category_colours(levels(df$group))
  }

  p <- ggplot(df, aes(x = group, y = value)) +
    geom_violin(aes(fill = group), trim = TRUE, alpha = 0.70, width = 1.05, position = position_identity()) +
    geom_boxplot(width = 0.07, outlier.alpha = 0.35, alpha = 0.90,
                 fill = "white", colour = "#222222", position = position_identity()) +
    annotate_top_right(sig_label(pval)) +
    scale_x_discrete(labels = labs_map) +
    scale_fill_manual(values = cols, guide = "none") +
    labs(
      x = NULL,
      y = y_label,
      title = metric_plain_label(metric),
      subtitle = paste0("Grouped by ", clean_plot_label(group), ". ", subtitle)
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.subtitle = element_text(lineheight = 0.95))

  save_plot_pair(p, stem, args$output_dir, 10.5, 6.3, formats)
  p + labs(title = metric_plain_label(metric))
}

p1 <- plot_violin_group("genome_gc", "phylum", "genome_GC_by_phylum", "Genome GC", collapse = TRUE)
p2 <- plot_violin_group("mean_GC3s", "phylum", "CDS_GC3s_by_phylum", "Mean CDS GC3s", collapse = TRUE)
p3 <- plot_violin_group("mean_tAI", "dikarya_status", "dikarya_non_dikarya_tAI_CDS", "Mean CDS tAI", collapse = FALSE)

# Diagnostic table for checking apparently low group values. A value below 0.5
# is biologically possible and does not by itself indicate an error; this table
# exposes sample size, range, mean and median and flags values outside [0, 1].
if ("mean_GC3s" %in% names(genomes)) {
  gc3_diag <- genomes %>%
    mutate(mean_GC3s = suppressWarnings(as.numeric(mean_GC3s)), phylum = clean_plot_label(phylum)) %>%
    filter(is.finite(mean_GC3s)) %>%
    group_by(phylum) %>%
    summarise(
      n_genomes = sum(is.finite(mean_GC3s)),
      group_mean_GC3s = mean(mean_GC3s, na.rm = TRUE),
      median_GC3s = median(mean_GC3s, na.rm = TRUE),
      min_GC3s = min(mean_GC3s, na.rm = TRUE),
      max_GC3s = max(mean_GC3s, na.rm = TRUE),
      any_outside_fraction_range = any(mean_GC3s < 0 | mean_GC3s > 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(mean_GC3s = group_mean_GC3s)
  readr::write_tsv(gc3_diag, file.path(args$output_dir, "CDS_GC3s_by_phylum_diagnostics.tsv"), na = "NA")
}

# GC and GC3s have the same 0-1 scale, so a split violin is a valid direct
# Dikarya/non-Dikarya comparison. Unrelated metrics are not mixed here.
if (all(c("genome_gc", "mean_GC3s", "dikarya_status") %in% names(genomes))) {
  split_gc <- bind_rows(
    genomes %>% transmute(metric = "Genome GC [%]", status = dikarya_status,
                          value = 100 * suppressWarnings(as.numeric(genome_gc))),
    genomes %>% transmute(metric = "mean_GC3s", status = dikarya_status,
                          value = 100 * suppressWarnings(as.numeric(mean_GC3s)))
  ) %>%
    filter(is.finite(value), !is.na(status)) %>%
    mutate(status = factor(status, levels = c("non-Dikarya", "Dikarya")),
           metric = factor(metric, levels = c("Genome GC [%]", "mean_GC3s")),
           split_side = split_side_from_level(status, "non-Dikarya"))
  counts <- split_gc %>% filter(metric == levels(metric)[1]) %>% count(status, name = "n")
  status_labels <- stats::setNames(paste0(counts$status, " (", counts$n, ")"), counts$status)
  p_split_gc <- ggplot(split_gc, aes(metric, value, fill = status, group = interaction(metric, status), split_side = split_side)) +
    geom_split_violin_project(trim = TRUE, alpha = 0.74, width = 0.93, scale = "width",
                              colour = "grey30", linewidth = 0.25) +
    geom_boxplot(aes(group = interaction(metric, status)), width = 0.08,
                 position = position_dodge(width = 0.13), outlier.shape = NA,
                 fill = "white", alpha = 0.88) +
    scale_fill_manual(values = binary_grey_values(levels(split_gc$status)), labels = status_labels, drop = FALSE) +
    scale_x_discrete(labels = c(
      "Genome GC [%]" = "Genome GC",
      "mean_GC3s" = "Mean CDS GC3s"
    )) +
    labs(x = NULL, y = "GC content [%]", fill = NULL,
         title = "Genome GC and mean CDS GC3s in Dikarya vs non-Dikarya",
         subtitle = "Left half: non-Dikarya; right half: Dikarya. Each observation is one genome.") +
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.text = element_text(size = 13.2),
      axis.title = element_text(size = 14.2),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 13),
      legend.text = element_text(size = 13),
      legend.key.width = grid::unit(1.0, "cm"),
      legend.spacing.x = grid::unit(0.45, "cm")
    )
  save_plot_pair(p_split_gc, "dikarya_non_dikarya_GC_GC3s_split_violin", args$output_dir, 10.5, 6.7, formats)
}

make_dikarya_panel <- function(metric, title, y_label, as_percent = FALSE) {
  if (!all(c(metric, "dikarya_status") %in% names(genomes))) return(NULL)
  d <- genomes %>%
    transmute(
      status = factor(dikarya_status, levels = c("non-Dikarya", "Dikarya")),
      value = suppressWarnings(as.numeric(.data[[metric]]))
    ) %>%
    filter(!is.na(status), is.finite(value))
  if (as_percent) d <- d %>% mutate(value = 100 * value)
  if (nrow(d) < 4L || n_distinct(d$status) < 2L) return(NULL)
  counts <- d %>% count(status, name = "n")
  labels <- stats::setNames(paste0(counts$status, " (", counts$n, ")"), counts$status)
  pval <- safe_wilcox_p(d$value, d$status)
  ggplot(d, aes(status, value, fill = status)) +
    geom_violin(trim = TRUE, alpha = 0.74, width = 0.92, scale = "width") +
    geom_boxplot(width = 0.10, outlier.shape = NA, fill = "white", alpha = 0.90) +
    annotate_top_right(sig_label(pval), size = 2.7) +
    scale_x_discrete(labels = labels) +
    scale_fill_manual(values = c("non-Dikarya" = "#B8B8B8", "Dikarya" = "#1F4E79"), guide = "none") +
    labs(x = NULL, y = y_label, title = title) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.x = element_text(angle = 18, hjust = 1, size = 13.2),
      axis.text.y = element_text(size = 13.4),
      axis.title.y = element_text(size = 14.4),
      plot.title = element_text(size = 15.8, face = "bold", hjust = 0.5),
      plot.margin = margin(6, 8, 8, 8)
    )
}

if (requireNamespace("patchwork", quietly = TRUE) && !is.null(p1) && !is.null(p2)) {
  top_grid <- patchwork::wrap_plots(list(p1, p2), ncol = 1)
  binary_panels <- Filter(Negate(is.null), list(
    make_dikarya_panel("mean_tAI", "Mean CDS tAI", "Mean CDS tAI"),
    make_dikarya_panel("fraction_lcr", "Genes with LCR", "Genes [%]", as_percent = TRUE),
    make_dikarya_panel("fraction_tm", "Genes with TM", "Genes [%]", as_percent = TRUE),
    make_dikarya_panel("fraction_signal_peptide", "Genes with signal peptide", "Genes [%]", as_percent = TRUE),
    make_dikarya_panel("fraction_pfam", "Genes with PFAM", "Genes [%]", as_percent = TRUE),
    make_dikarya_panel("fraction_with_go", "Genes with GO annotation", "Genes [%]", as_percent = TRUE)
  ))
  if (length(binary_panels) > 0L) {
    bottom_grid <- patchwork::wrap_plots(binary_panels, ncol = 2)
    grid <- patchwork::wrap_plots(list(top_grid, bottom_grid), ncol = 1, heights = c(2, 3)) +
      patchwork::plot_annotation(
        title = "Genome-level GC metrics and Dikarya/non-Dikarya feature summaries"
      )
  } else {
    grid <- top_grid
  }
  save_plot_pair(grid, "genome_metric_violins_grid", args$output_dir, 15.0, 26.0, formats)
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
      annotate_top_right(ann, size = 3.6) +
      scale_colour_manual(values = c("non-Dikarya" = "#B8B8B8", "Dikarya" = "#1F4E79"), drop = FALSE) +
      labs(
        x = "Genome GC",
        y = "Mean CDS GC3s",
        colour = NULL,
        title = "Genome GC vs mean CDS GC3s"
      ) +
      theme_minimal(base_size = 12)
    save_plot_pair(p, "genome_GC_vs_CDS_GC3s", args$output_dir, 8, 5.7, formats)
  }
}

# Missing legacy information retained as one non-redundant panel: mean CDS tAI
# is related to mean CDS GC and GC3s, with Dikarya status shown consistently.
gc_predictors <- intersect(c("mean_GC", "mean_GC3s"), names(genomes))
if ("mean_tAI" %in% names(genomes) && length(gc_predictors) > 0L) {
  tai_gc <- purrr::map_dfr(gc_predictors, function(predictor) {
    genomes %>%
      transmute(
        sample = if ("sample" %in% names(genomes)) as.character(sample) else as.character(row_number()),
        dikarya_status,
        predictor = predictor,
        gc_percent = 100 * suppressWarnings(as.numeric(.data[[predictor]])),
        mean_tAI = suppressWarnings(as.numeric(mean_tAI))
      )
  }) %>%
    filter(is.finite(gc_percent), is.finite(mean_tAI), !is.na(dikarya_status)) %>%
    mutate(
      dikarya_status = factor(dikarya_status, levels = c("non-Dikarya", "Dikarya")),
      predictor = factor(
        predictor,
        levels = c("mean_GC", "mean_GC3s"),
        labels = c("Mean CDS GC [%]", "Mean CDS GC3s [%]")
      )
    )
  if (nrow(tai_gc) >= 6L) {
    regression_summary <- tai_gc %>%
      group_by(predictor, dikarya_status) %>%
      group_modify(~ {
        fit <- if (nrow(.x) >= 3L && n_distinct(.x$gc_percent) >= 2L) stats::lm(mean_tAI ~ gc_percent, data = .x) else NULL
        if (is.null(fit)) return(tibble(n_genomes = nrow(.x), intercept = NA_real_, slope = NA_real_, r_squared = NA_real_, p_value = NA_real_))
        sm <- summary(fit)
        tibble(
          n_genomes = nrow(.x), intercept = unname(coef(fit)[[1]]),
          slope = unname(coef(fit)[[2]]),
          r_squared = sm$r.squared, p_value = sm$coefficients[2, 4]
        )
      }) %>%
      ungroup()

    regression_labels <- regression_summary %>%
      mutate(
        label = sprintf(
          "y = %.3f %+.4fx; R² = %.3f; p = %s",
          intercept, slope, r_squared,
          format.pval(p_value, digits = 2, eps = 0.001)
        ),
        vjust = if_else(as.character(dikarya_status) == "non-Dikarya", 1.1, 2.5)
      )

    p_tai_gc <- ggplot(tai_gc, aes(gc_percent, mean_tAI, colour = dikarya_status)) +
      geom_point(size = 2.2, alpha = 0.82) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.78, alpha = 0.09) +
      geom_text(
        data = regression_labels,
        aes(x = -Inf, y = Inf, label = label, colour = dikarya_status, vjust = vjust),
        inherit.aes = FALSE, hjust = -0.03, size = 3.0, show.legend = FALSE
      ) +
      facet_wrap(vars(predictor), nrow = 1, scales = "free_x") +
      scale_colour_manual(
        values = c("non-Dikarya" = "#B8B8B8", "Dikarya" = "#1F4E79"),
        drop = FALSE
      ) +
      labs(
        x = NULL, y = "Mean CDS tAI", colour = NULL,
        title = "Mean CDS tAI versus CDS GC composition",
        subtitle = "Lines and confidence bands are fitted separately for Dikarya and non-Dikarya genomes."
      ) +
      theme_minimal(base_size = 15) +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text = element_text(size = 13.2),
        axis.title = element_text(size = 14.2),
        strip.text = element_text(face = "bold", size = 14.2),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.key.width = grid::unit(1.0, "cm"),
        legend.spacing.x = grid::unit(0.45, "cm")
      )
    save_plot_pair(p_tai_gc, "dikarya_tai_vs_cds_gc_gc3", args$output_dir, 12.0, 6.7, formats)

    readr::write_tsv(
      regression_summary,
      file.path(args$output_dir, "dikarya_tai_vs_cds_gc_gc3_regression.tsv"),
      na = "NA"
    )
  }
}

plot_tai_predictors <- function(metrics, labels, stem, title, log10_x = FALSE) {
  available <- intersect(metrics, names(genomes))
  if (length(available) == 0L) {
    add_missing(stem, paste("missing all predictors:", paste(metrics, collapse = ", ")))
    return(NULL)
  }

  long <- genomes %>%
    select(mean_tAI, any_of(available)) %>%
    pivot_longer(cols = any_of(available), names_to = "predictor", values_to = "value") %>%
    mutate(
      mean_tAI = suppressWarnings(as.numeric(mean_tAI)),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    filter(is.finite(mean_tAI), is.finite(value), !log10_x | value > 0) %>%
    mutate(
      x_plot = if (log10_x) log10(value) else value,
      predictor_label = factor(predictor, levels = available, labels = unname(labels[available]))
    )

  if (nrow(long) == 0L) {
    add_missing(stem, "no finite predictor/tAI pairs")
    return(NULL)
  }

  p <- ggplot(long, aes(x_plot, mean_tAI)) +
    geom_point(size = 1.8, alpha = 0.68, colour = "#1F4E79") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.7, colour = "#4D4D4D") +
    facet_wrap(vars(predictor_label), scales = "free_x", ncol = 2) +
    labs(
      x = if (log10_x) "Predictor (log10 scale)" else "Predictor",
      y = "Mean CDS tAI",
      title = title,
      subtitle = "Each observation is one genome; lines show separate descriptive linear fits."
    ) +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

  save_plot_pair(p, stem, args$output_dir, 10.5, 7.5, formats)
  long
}

assembly_data <- plot_tai_predictors(
  c("n50_bp", "l50", "n_contigs"),
  c(n50_bp = "Assembly N50 [bp]", l50 = "Assembly L50", n_contigs = "Number of contigs"),
  "tai_vs_assembly_quality",
  "Mean CDS tAI versus assembly quality statistics",
  log10_x = TRUE
)
scale_data <- plot_tai_predictors(
  c("genome_size_bp", "n_genes", "n_proteins", "proteome_length_aa"),
  c(
    genome_size_bp = "Genome size [bp]", n_genes = "Number of CDS",
    n_proteins = "Number of proteins", proteome_length_aa = "Proteome length [aa]"
  ),
  "tai_vs_genome_scale",
  "Mean CDS tAI versus genome and proteome scale",
  log10_x = TRUE
)
gc_data <- plot_tai_predictors(
  c("genome_gc", "mean_GC3s"),
  c(genome_gc = "Genome GC", mean_GC3s = "Mean CDS GC3s"),
  "tai_vs_genome_gc",
  "Mean CDS tAI versus genomic GC composition",
  log10_x = FALSE
)

predictor_data <- bind_rows(assembly_data, scale_data, gc_data)
if (nrow(predictor_data) > 0L) {
  predictor_correlations <- predictor_data %>%
    group_by(predictor) %>%
    summarise(
      n_genomes = n(),
      spearman_rho = if (n() >= 3L && n_distinct(value) >= 2L) {
        suppressWarnings(stats::cor(value, mean_tAI, method = "spearman"))
      } else NA_real_,
      p_value = if (n() >= 3L && n_distinct(value) >= 2L) {
        suppressWarnings(stats::cor.test(value, mean_tAI, method = "spearman", exact = FALSE)$p.value)
      } else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))
  readr::write_tsv(
    predictor_correlations,
    file.path(args$output_dir, "tai_genome_feature_correlations.tsv"),
    na = "NA"
  )
}
if (!"busco_score" %in% names(genomes)) {
  add_missing("tai_vs_busco", "BUSCO score is not present in genome_summary.tsv")
}

if (nrow(missing) == 0L) missing <- tibble(plot = "none", reason = "all possible genome GC/tAI plots generated")
readr::write_tsv(missing, file.path(args$output_dir, "XXXXX_missing_plots.tsv"))
message("Done: ", args$output_dir)
