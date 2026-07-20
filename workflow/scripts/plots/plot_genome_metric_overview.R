#!/usr/bin/env Rscript
# Render genome-level metric summaries by phylum and lifestyle.
# This script is intentionally based on genome_metric_summary.tsv, not raw gene rows.
# It collapses rare categories into Other and avoids reusing column name `n` before
# weighted summaries are calculated.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

source("workflow/scripts/lib/plot_style_helpers.R")

safe_read <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

save_plot <- function(plot, path, width = 10, height = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = 300, bg = "white", limitsize = FALSE)
}

sig_label <- function(p) {
  p_label(p)[[1]]
}

collapse_summary <- function(data, min_n = 5L) {
  required <- c("group_variable", "group_value", "metric", "n", "mean", "median", "q1", "q3")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) stop("Metric summary lacks columns: ", paste(missing, collapse = ", "))

  data %>%
    mutate(
      n = suppressWarnings(as.numeric(n)),
      mean = suppressWarnings(as.numeric(mean)),
      median = suppressWarnings(as.numeric(median)),
      q1 = suppressWarnings(as.numeric(q1)),
      q3 = suppressWarnings(as.numeric(q3)),
      group_value_raw = clean_label(group_value),
      group_value_collapsed = if_else(n >= min_n, group_value_raw, "Other"),
      .weight = if_else(is.finite(n) & n > 0, n, 0)
    ) %>%
    group_by(group_variable, group_value_collapsed, metric) %>%
    summarise(
      n = sum(.weight, na.rm = TRUE),
      mean = if (sum(.weight, na.rm = TRUE) > 0) stats::weighted.mean(mean, w = .weight, na.rm = TRUE) else NA_real_,
      median = if (sum(.weight, na.rm = TRUE) > 0) stats::weighted.mean(median, w = .weight, na.rm = TRUE) else NA_real_,
      q1 = if (sum(.weight, na.rm = TRUE) > 0) stats::weighted.mean(q1, w = .weight, na.rm = TRUE) else NA_real_,
      q3 = if (sum(.weight, na.rm = TRUE) > 0) stats::weighted.mean(q3, w = .weight, na.rm = TRUE) else NA_real_,
      other_members = paste(sort(unique(group_value_raw[group_value_collapsed == "Other"])), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      group_label = clean_label(group_value_collapsed),
      group_label_n = paste0(group_label, " (", as.integer(n), ")"),
      metric_label = metric_plain_label(metric)
    )
}

get_group_p <- function(effects, group_name, metric_name) {
  if (!all(c("group_variable", "metric", "p_value") %in% names(effects))) return(NA_real_)
  hit <- effects %>%
    filter(.data$group_variable == group_name, .data$metric == metric_name) %>%
    pull(.data$p_value)
  hit <- suppressWarnings(as.numeric(hit))
  hit <- hit[is.finite(hit)]
  if (length(hit) == 0L) NA_real_ else hit[[1]]
}

write_metric_range_diagnostics <- function(d, group_name, output_dir) {
  diag <- d %>%
    group_by(metric, metric_label) %>%
    summarise(
      min_q1 = min(q1, median, na.rm = TRUE),
      max_q3 = max(q3, median, na.rm = TRUE),
      range = max_q3 - min_q1,
      min_median = min(median, na.rm = TRUE),
      max_median = max(median, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(max_q3))
  readr::write_tsv(diag, file.path(output_dir, paste0("genome_metric_range_diagnostics_by_", group_name, ".tsv")), na = "NA")
}

make_group_plot <- function(metrics, effects, group_name, title, output_dir, min_n = 5L) {
  d <- metrics %>%
    filter(.data$group_variable == group_name, !stringr::str_detect(.data$metric, "^median_")) %>%
    collapse_summary(min_n = min_n) %>%
    filter(is.finite(.data$mean))

  if (nrow(d) == 0L) stop("No finite data for group variable: ", group_name)
  write_metric_range_diagnostics(d, group_name, output_dir)

  group_levels <- d %>%
    distinct(group_value_collapsed, group_label_n) %>%
    mutate(group_value_collapsed = factor(group_value_collapsed, levels = order_project_groups(group_value_collapsed, group_name))) %>%
    arrange(group_value_collapsed, group_label_n) %>%
    pull(group_label_n)

  d <- d %>%
    mutate(
      group_label_n = factor(group_label_n, levels = unique(group_levels)),
      metric_label = factor(
        metric_label,
        levels = unique(metric_plain_label(metrics$metric[!stringr::str_detect(metrics$metric, "^median_")]))
      )
    )

  other <- d %>%
    filter(group_value_collapsed == "Other", nzchar(other_members)) %>%
    pull(other_members) %>%
    unique()
  subtitle <- if (length(other) > 0L) {
    paste0("Other contains categories with n < ", min_n, ": ", paste(clean_label(unlist(strsplit(other, "; "))), collapse = ", "))
  } else {
    NULL
  }

  pvals <- d %>%
    distinct(metric) %>%
    rowwise() %>%
    mutate(p_label = sig_label(get_group_p(effects, group_name, metric))) %>%
    ungroup()

  global_means <- d %>%
    group_by(metric) %>%
    summarise(global_mean = stats::weighted.mean(mean, w = n, na.rm = TRUE), .groups = "drop")

  priority <- c("mean_tAI", "mean_CAI", "genome_gc", "mean_GC3s")
  metric_order <- c(intersect(priority, unique(d$metric)), sort(setdiff(unique(d$metric), priority)))
  d <- d %>%
    mutate(
      metric = factor(metric, levels = metric_order),
      scale_group = case_when(
        as.character(metric) %in% c("mean_tAI", "mean_CAI") ~ "codon_adaptation",
        as.character(metric) %in% c("genome_gc", "mean_GC3s") ~ "gc_content",
        TRUE ~ as.character(metric)
      )
    )

  range_table <- d %>%
    group_by(scale_group) %>%
    summarise(
      y_min = min(c(mean, q1, q3), na.rm = TRUE),
      y_max = max(c(mean, q1, q3), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      span = pmax(y_max - pmin(0, y_min), 1e-6),
      y_min_plot = pmin(0, y_min - 0.04 * span),
      y_max_plot = y_max + 0.08 * span
    )

  colour_lookup <- d %>%
    distinct(group_label, group_label_n)
  category_colours <- project_category_colours(colour_lookup$group_label)
  colour_lookup <- colour_lookup %>%
    mutate(colour = unname(category_colours[as.character(group_label)]))
  colour_map <- stats::setNames(colour_lookup$colour, colour_lookup$group_label_n)

  make_metric_panel <- function(metric_name) {
    dm <- d %>% filter(as.character(metric) == metric_name)
    limits <- range_table %>% filter(scale_group == first(dm$scale_group))
    annotation <- pvals %>%
      filter(metric == metric_name) %>%
      left_join(global_means %>% filter(metric == metric_name), by = "metric")

    ggplot(dm, aes(x = group_label_n, y = mean, fill = group_label_n)) +
      geom_col(width = 0.60, alpha = 0.88, colour = "grey30", linewidth = 0.15) +
      geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.18, linewidth = 0.45) +
      geom_hline(
        yintercept = annotation$global_mean[[1]],
        linetype = "dashed", linewidth = 0.70, colour = "#1F4E79"
      ) +
      scale_x_discrete(labels = function(x) stringr::str_wrap(x, 18)) +
      scale_y_continuous(limits = c(limits$y_min_plot[[1]], limits$y_max_plot[[1]]), expand = expansion(mult = 0)) +
      scale_fill_manual(values = colour_map, guide = "none") +
      labs(
        x = NULL, y = NULL,
        title = metric_plain_label(metric_name),
        tag = annotation$p_label[[1]]
      ) +
      theme_minimal(base_size = 12.5) +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 10.2),
        axis.text.y = element_text(size = 10.5),
        plot.title = element_text(size = 12.8, face = "bold", hjust = 0),
        plot.tag = element_text(size = 10.5, hjust = 1, face = "plain"),
        plot.tag.position = c(0.98, 0.99),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(9, 12, 12, 12)
      )
  }

  panels <- lapply(metric_order, make_metric_panel)
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for the paired genome-metric layout.")
  }
  patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = title,
      subtitle = if (is.null(subtitle)) NULL else wrap_text(subtitle, 105),
      caption = "Bars show group means; error bars show the IQR of genome-level values; the dashed line is the weighted mean across displayed groups."
    ) &
    theme(
      plot.title = element_text(size = 17, face = "bold"),
      plot.subtitle = element_text(size = 11, lineheight = 1.15),
      plot.caption = element_text(size = 10.2),
      plot.margin = margin(14, 16, 14, 16)
    )
}

if (exists("snakemake")) {
  metric_path <- snakemake@input[["metric_summary"]]
  effect_path <- snakemake@input[["effect_summary"]]
  metrics <- safe_read(metric_path)
  effects <- safe_read(effect_path)

  outdir <- dirname(snakemake@output[["phylum_png"]])
  phylum_plot <- make_group_plot(metrics, effects, "phylum", "Genome metrics by phylum", outdir, min_n = 5L)
  lifestyle_plot <- make_group_plot(metrics, effects, "lifestyle", "Genome metrics by lifestyle", outdir, min_n = 5L)

  message("Saving genome-metric plots by phylum")
  save_plot(phylum_plot, snakemake@output[["phylum_png"]], width = 12.2, height = 15.0)
  save_plot(phylum_plot, snakemake@output[["phylum_pdf"]], width = 12.2, height = 15.0)
  message("Saving genome-metric plots by lifestyle")
  save_plot(lifestyle_plot, snakemake@output[["lifestyle_png"]], width = 12.4, height = 15.2)
  save_plot(lifestyle_plot, snakemake@output[["lifestyle_pdf"]], width = 12.4, height = 15.2)
} else {
  stop("This script is intended to be run by Snakemake.")
}
