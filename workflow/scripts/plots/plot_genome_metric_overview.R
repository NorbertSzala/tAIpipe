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

clean_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- stringr::str_squish(x)
  stringr::str_to_sentence(x)
}

wrap_text <- function(x, width = 95) stringr::str_wrap(x, width = width)

safe_read <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

save_plot <- function(plot, path, width = 10, height = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
}

sig_label <- function(p) {
  if (!is.finite(p)) return("p = NA")
  paste0("p = ", format.pval(p, digits = 2, eps = 1e-4), ifelse(p < 0.05, "\nsignificant", "\nnot significant"))
}

annotate_top_right <- function(label) {
  ggplot2::annotate(
    "text", x = Inf, y = Inf, label = label,
    hjust = 1.05, vjust = 1.15, size = 3.0
  )
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
      group_value_raw = as.character(group_value),
      group_value_collapsed = if_else(n > min_n, group_value_raw, "Other"),
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
      metric_label = clean_label(metric)
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

make_group_plot <- function(metrics, effects, group_name, title, min_n = 5L) {
  d <- metrics %>%
    filter(.data$group_variable == group_name) %>%
    collapse_summary(min_n = min_n) %>%
    filter(is.finite(.data$median))

  if (nrow(d) == 0L) stop("No finite data for group variable: ", group_name)

  # Order categories by median across all displayed metrics.
  order_tbl <- d %>%
    group_by(group_value_collapsed, group_label_n) %>%
    summarise(order_value = median(median, na.rm = TRUE), .groups = "drop") %>%
    arrange(order_value)

  d <- d %>%
    mutate(group_label_n = factor(group_label_n, levels = order_tbl$group_label_n))

  other <- d %>%
    filter(group_value_collapsed == "Other", nzchar(other_members)) %>%
    pull(other_members) %>%
    unique()
  subtitle <- if (length(other) > 0L) {
    wrap_text(paste0("Other contains categories with n <= ", min_n, ": ", paste(clean_label(unlist(strsplit(other, "; "))), collapse = ", ")), 120)
  } else {
    NULL
  }

  # One global y-scale for all metrics in this figure. These metrics are fractions/indexes on similar scales.
  y_min <- min(d$q1, d$median, na.rm = TRUE)
  y_max <- max(d$q3, d$median, na.rm = TRUE)
  pad <- 0.06 * (y_max - y_min)
  if (!is.finite(pad) || pad == 0) pad <- 0.05

  pvals <- d %>%
    distinct(metric, metric_label) %>%
    rowwise() %>%
    mutate(p_label = sig_label(get_group_p(effects, group_name, metric))) %>%
    ungroup()

  d <- d %>% left_join(pvals, by = c("metric", "metric_label"))

  ggplot(d, aes(x = group_label_n, y = median, fill = group_label_n)) +
    geom_col(width = 0.72, alpha = 0.88, colour = "grey30", linewidth = 0.15) +
    geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.22, linewidth = 0.45) +
    geom_text(
      data = pvals,
      aes(x = Inf, y = Inf, label = p_label),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = 1.15, size = 2.7
    ) +
    facet_wrap(vars(metric_label), scales = "fixed", ncol = 1) +
    coord_cartesian(ylim = c(y_min - pad, y_max + 1.8 * pad), clip = "off") +
    labs(
      x = NULL,
      y = "Median genome-level value (IQR)",
      title = title,
      subtitle = subtitle,
      caption = "Bars show median of genome-level values; error bars show IQR. Category labels include number of genomes."
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      plot.subtitle = element_text(lineheight = 0.95),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(7, 18, 7, 7)
    )
}

if (exists("snakemake")) {
  metric_path <- snakemake@input[["metric_summary"]]
  effect_path <- snakemake@input[["effect_summary"]]
  metrics <- safe_read(metric_path)
  effects <- safe_read(effect_path)

  phylum_plot <- make_group_plot(metrics, effects, "phylum", "Genome metrics by phylum", min_n = 5L)
  lifestyle_plot <- make_group_plot(metrics, effects, "lifestyle", "Genome metrics by lifestyle", min_n = 5L)

  save_plot(phylum_plot, snakemake@output[["phylum_png"]], width = 8.5, height = 9.0)
  save_plot(phylum_plot, snakemake@output[["phylum_pdf"]], width = 8.5, height = 9.0)
  save_plot(lifestyle_plot, snakemake@output[["lifestyle_png"]], width = 9.5, height = 9.5)
  save_plot(lifestyle_plot, snakemake@output[["lifestyle_pdf"]], width = 9.5, height = 9.5)
} else {
  stop("This script is intended to be run by Snakemake.")
}
