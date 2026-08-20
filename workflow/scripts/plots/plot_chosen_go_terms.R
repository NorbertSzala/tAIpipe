#!/usr/bin/env Rscript
# Data-driven plots for the selected GO terms.
#
# Interpretation notes for reviewers and future users:
# - The enrichment effect is the Cochran-Mantel-Haenszel (CMH) common odds ratio
#   across genomes. Each genome forms a stratum with a 2x2 table: gene is in the
#   selected tAI tail vs outside the tail, and gene has the GO term vs does not
#   have the GO term. The CMH common odds ratio estimates one genome-stratified
#   enrichment effect while reducing confounding by genome size, annotation depth
#   and lineage-specific gene counts. log2(CMH common OR) is plotted so that 0
#   means no enrichment, positive values mean enrichment in the selected tail,
#   and negative values mean depletion from the selected tail.
# - -log10(FDR) is a monotonic transformation of the multiple-testing-adjusted
#   q-value. Larger values mean stronger statistical support after FDR correction.
#   It is used only for visual point size/significance ranking, not as an effect
#   size. q=0 values from numerical underflow are clipped to the smallest positive
#   double before taking -log10.
# - tAI_z is within-genome standardized tAI. It is preferred for cross-genome
#   distribution plots because raw tAI scales can differ between genomes.

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
    chosen_go_terms = "results/statistics/chosen_GOterms.tsv",
    gene_table = "results/tables/gene_features.tsv",
    output_dir = "results/plots/go_chosen_terms",
    formats = "png,pdf",
    max_terms = "48",
    tail_fraction = "0.10",
    max_tail_genes = "0",
    min_genes_with_go_per_sample = "10",
    manifest_output = "results/plots/go_chosen_terms/plot_manifest.tsv"
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
tail_fraction <- as.numeric(args$tail_fraction)
max_tail_genes <- as.integer(args$max_tail_genes)
min_genes_with_go_per_sample <- as.integer(args$min_genes_with_go_per_sample)
if (!is.finite(tail_fraction) || tail_fraction <= 0 || tail_fraction >= 0.5) stop("tail_fraction must be in (0, 0.5)")
if (!is.finite(max_tail_genes) || max_tail_genes < 0L) stop("max_tail_genes must be >= 0")
if (!is.finite(min_genes_with_go_per_sample) || min_genes_with_go_per_sample < 1L) stop("min_genes_with_go_per_sample must be >= 1")
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(args$output_dir, c(
  "chosen_go_terms_enrichment_effect_barplot.png", "chosen_go_terms_enrichment_effect_barplot.pdf",
  "chosen_go_terms_support_effect_bubble.png", "chosen_go_terms_support_effect_bubble.pdf"
)), force = TRUE)

save_plot_pair_local <- function(p, stem, width = 9, height = 7) {
  p <- prepare_plot_for_export(p)
  for (fmt in formats) {
    ggplot2::ggsave(
      file.path(args$output_dir, paste0(stem, ".", fmt)),
      p, width = width, height = height, units = "in", dpi = 300, bg = "white", limitsize = FALSE
    )
  }
}

write_plot_manifest <- function() {
  manifest_files <- list.files(
    args$output_dir,
    pattern = "\\.(png|pdf|tsv)$",
    full.names = TRUE
  )
  manifest <- tibble::tibble(
    file = basename(manifest_files),
    size_bytes = file.info(manifest_files)$size
  ) %>% arrange(file)
  readr::write_tsv(manifest, args$manifest_output)
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  stringr::str_squish(x)
}

clean_namespace <- function(x) {
  clean_label(x)
}

clean_tail <- function(x) {
  x <- clean_text(x)
  dplyr::case_when(
    x %in% c("high_tAI", "High tai", "High tAI") ~ "High tAI",
    x %in% c("low_tAI", "Low tai", "Low tAI") ~ "Low tAI",
    TRUE ~ clean_label(x)
  )
}

first_existing <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) NA_character_ else hit[[1]]
}

as_num <- function(x) suppressWarnings(as.numeric(x))

pick_numeric_effect <- function(df) {
  log2_cols <- c("log2_common_odds_ratio", "common_log2_odds_ratio", "cmh_log2_or", "log2_odds_ratio", "log2_or", "log2OR")
  log_cols <- c("common_log_odds_ratio", "cmh_log_or", "log_odds_ratio", "log_or")
  or_cols <- c("common_odds_ratio", "cmh_odds_ratio", "common_or", "odds_ratio", "OR", "estimate")

  for (col in intersect(log2_cols, names(df))) {
    val <- as_num(df[[col]])
    if (any(is.finite(val))) return(list(column = col, type = "log2", value = val))
  }
  for (col in intersect(log_cols, names(df))) {
    val <- as_num(df[[col]]) / log(2)
    if (any(is.finite(val))) return(list(column = col, type = "natural_log", value = val))
  }
  for (col in intersect(or_cols, names(df))) {
    raw <- as_num(df[[col]])
    val <- ifelse(is.finite(raw) & raw > 0, log2(raw), NA_real_)
    if (any(is.finite(val))) return(list(column = col, type = "odds_ratio", value = val))
  }
  list(column = NA_character_, type = NA_character_, value = rep(NA_real_, nrow(df)))
}

extract_go_ids <- function(x) {
  ids <- stringr::str_extract_all(toupper(as.character(x)), "GO:[0-9]{7}")
  ids <- unlist(ids, use.names = FALSE)
  unique(ids[!is.na(ids) & nzchar(ids)])
}

missing <- tibble(plot = character(), reason = character())
add_missing <- function(plot, reason) {
  missing <<- bind_rows(missing, tibble(plot = plot, reason = reason))
  message("XXXXX skipped ", plot, ": ", reason)
}

if (!file.exists(args$chosen_go_terms)) stop("chosen_GOterms.tsv does not exist: ", args$chosen_go_terms)
chosen <- readr::read_tsv(args$chosen_go_terms, show_col_types = FALSE, progress = FALSE)

if (nrow(chosen) == 0L) {
  add_missing("all_go_plots", "chosen_GOterms.tsv has zero rows")
  readr::write_tsv(missing, file.path(args$output_dir, "XXXXX_missing_go_plots.tsv"))
  write_plot_manifest()
  quit(save = "no", status = 0)
}

chosen <- chosen %>%
  mutate(across(where(is.character), clean_text)) %>%
  filter(if ("go_id" %in% names(.)) go_id != "XXXXX" else TRUE)

if (nrow(chosen) == 0L) {
  add_missing("all_go_plots", "chosen_GOterms.tsv contains only XXXXX placeholder rows")
  readr::write_tsv(missing, file.path(args$output_dir, "XXXXX_missing_go_plots.tsv"))
  write_plot_manifest()
  quit(save = "no", status = 0)
}

go_col <- first_existing(chosen, c("go_id", "go_term", "go", "term", "go_terms"))
name_col <- first_existing(chosen, c("go_name", "term_name", "name", "description"))
namespace_col <- first_existing(chosen, c("go_namespace", "namespace", "ontology"))
tail_col <- first_existing(chosen, c("tail", "tai_tail", "direction", "contrast", "comparison", "group"))
q_col <- first_existing(chosen, c("q_value", "q", "fdr", "padj", "adjusted_p_value", "p_adjust", "p_adj"))
p_col <- first_existing(chosen, c("p_value", "p", "pvalue"))
support_col <- first_existing(chosen, c("n_informative_genomes", "n_genomes_with_term", "n_samples_with_go_in_gene_features"))

if (is.na(go_col)) stop("chosen_GOterms.tsv must contain a GO ID column, e.g. go_id")

effect <- pick_numeric_effect(chosen)

plot_df <- chosen %>%
  mutate(
    go_id_plot = as.character(.data[[go_col]]),
    go_name_plot = if (!is.na(name_col)) clean_text(.data[[name_col]]) else go_id_plot,
    namespace_plot = if (!is.na(namespace_col)) clean_namespace(.data[[namespace_col]]) else "GO",
    tail_plot = if (!is.na(tail_col)) clean_tail(.data[[tail_col]]) else "Selected GO terms",
    q_plot = if (!is.na(q_col)) as_num(.data[[q_col]]) else if (!is.na(p_col)) as_num(.data[[p_col]]) else NA_real_,
    p_plot = if (!is.na(p_col)) as_num(.data[[p_col]]) else NA_real_,
    support_plot = if (!is.na(support_col)) as_num(.data[[support_col]]) else NA_real_,
    log2_or_plot = effect$value,
    q_for_rank = ifelse(is.finite(q_plot), pmax(q_plot, .Machine$double.xmin), NA_real_),
    rank_score = ifelse(is.finite(q_for_rank), -log10(q_for_rank), NA_real_),
    rank_score = pmin(rank_score, 320),
    term_label_plain = paste0(go_id_plot, " - ", go_name_plot)
  ) %>%
  arrange(q_plot, desc(abs(log2_or_plot))) %>%
  slice_head(n = max_terms)

readr::write_tsv(
  tibble(
    detected_go_column = go_col,
    detected_name_column = name_col,
    detected_namespace_column = namespace_col,
    detected_tail_column = tail_col,
    detected_q_column = q_col,
    detected_p_column = p_col,
    detected_support_column = support_col,
    detected_effect_column = effect$column,
    detected_effect_type = effect$type,
    finite_effect_rows = sum(is.finite(plot_df$log2_or_plot))
  ),
  file.path(args$output_dir, "chosen_go_terms_detected_columns.tsv")
)

readr::write_tsv(
  plot_df %>% select(-q_for_rank) %>% mutate(across(where(is.character), clean_text)),
  file.path(args$output_dir, "chosen_go_terms_used_for_plots.tsv"),
  na = "NA"
)


tail_shapes <- c("High tAI" = 21, "Low tAI" = 24)
namespace_colours <- c(
  "biological process" = "#252525",
  "molecular function" = "#777777",
  "cellular component" = "#D0D0D0",
  "cellular components" = "#D0D0D0",
  "Biological process" = "#252525",
  "Molecular function" = "#777777",
  "Cellular component" = "#D0D0D0",
  "Cellular components" = "#D0D0D0"
)
tai_namespace_colours <- c(
  "biological process" = "#08306B",
  "molecular function" = "#2171B5",
  "cellular component" = "#6BAED6",
  "cellular components" = "#6BAED6",
  "Biological process" = "#08306B",
  "Molecular function" = "#2171B5",
  "Cellular component" = "#6BAED6",
  "Cellular components" = "#6BAED6"
)
namespace_facet_labeller <- labeller(
  namespace_plot = function(x) stringr::str_wrap(x, width = 12)
)

# Add tail-specific descriptive quantities from gene_features.tsv. These are
# descriptive only; the CMH odds ratio remains the enrichment effect estimate.
tail_metrics <- tibble()
if (file.exists(args$gene_table)) {
  genes_for_size <- readr::read_tsv(args$gene_table, show_col_types = FALSE, progress = FALSE)
  if (all(c("sample", "gene_id", "go_terms", "tAI") %in% names(genes_for_size))) {
    if (!"tAI_z" %in% names(genes_for_size)) genes_for_size$tAI_z <- NA_real_
    selected_ids <- unique(plot_df$go_id_plot)
    ranked_for_size <- genes_for_size %>%
      transmute(
        sample = as.character(sample), gene_id = as.character(gene_id),
        go_terms = as.character(go_terms),
        tAI = suppressWarnings(as.numeric(tAI)),
        tAI_z = suppressWarnings(as.numeric(tAI_z)),
        go_id_list = stringr::str_extract_all(toupper(as.character(go_terms)), "GO:[0-9]{7}")
      ) %>%
      filter(
        is.finite(tAI), !is.na(sample), nzchar(sample),
        lengths(go_id_list) > 0L
      ) %>%
      group_by(sample) %>%
      filter(n() >= min_genes_with_go_per_sample) %>%
      arrange(desc(tAI), gene_id, .by_group = TRUE) %>%
      mutate(rank_high = row_number()) %>%
      arrange(tAI, gene_id, .by_group = TRUE) %>%
      mutate(
        rank_low = row_number(),
        n_ranked = n(),
        tail_n_fraction = ceiling(tail_fraction * n_ranked),
        tail_n_uncapped = if (max_tail_genes > 0L) pmin(tail_n_fraction, max_tail_genes) else tail_n_fraction,
        tail_n = pmin(tail_n_uncapped, floor(n_ranked / 2)),
        high_tail = rank_high <= tail_n,
        low_tail = rank_low <= tail_n
      ) %>%
      ungroup()

    tail_metrics <- bind_rows(
      ranked_for_size %>% filter(high_tail) %>% mutate(tail_plot = "High tAI"),
      ranked_for_size %>% filter(low_tail) %>% mutate(tail_plot = "Low tAI")
    ) %>%
      tidyr::unnest_longer(go_id_list, values_to = "go_id_plot") %>%
      filter(go_id_plot %in% selected_ids) %>%
      distinct(sample, gene_id, go_id_plot, tail_plot, .keep_all = TRUE) %>%
      group_by(go_id_plot, tail_plot) %>%
      summarise(
        mean_tail_tAI = mean(tAI, na.rm = TRUE),
        mean_tail_tAI_z = mean(tAI_z, na.rm = TRUE),
        n_tail_genes = n_distinct(paste(sample, gene_id, sep = "\t")),
        n_tail_genomes = n_distinct(sample),
        .groups = "drop"
      )
    plot_df <- plot_df %>% left_join(tail_metrics, by = c("go_id_plot", "tail_plot"))
    readr::write_tsv(tail_metrics, file.path(args$output_dir, "chosen_go_terms_tail_descriptive_metrics.tsv"), na = "NA")
  }
}

prepare_tai_bar_overlay <- function(data, bar_col) {
  bar_values <- suppressWarnings(as.numeric(data[[bar_col]]))
  effect_limit <- max(abs(data$log2_or_plot[is.finite(data$log2_or_plot)]), na.rm = TRUE)
  bar_limit <- max(abs(bar_values[is.finite(bar_values)]), na.rm = TRUE)
  if (!is.finite(effect_limit) || effect_limit <= 0) effect_limit <- 1
  if (!is.finite(bar_limit) || bar_limit <= 0) bar_limit <- 1
  scale_factor <- effect_limit / bar_limit
  list(
    data = data %>% mutate(bar_value = bar_values, bar_scaled = bar_value * scale_factor),
    scale_factor = scale_factor,
    breaks = pretty(range(bar_values[is.finite(bar_values)], na.rm = TRUE), n = 4)
  )
}

tai_bar_spec <- function(size_col) {
  if (identical(size_col, "mean_tail_tAI") && "mean_tail_tAI_z" %in% names(plot_df)) {
    list(column = "mean_tail_tAI_z", title = "Mean tAI z-score")
  } else if ("mean_tail_tAI" %in% names(plot_df)) {
    list(column = "mean_tail_tAI", title = "Mean raw tAI")
  } else {
    NULL
  }
}

make_enrichment_dotplot <- function(size_col, size_title, stem, size_range = c(3.2, 9.5), abs_size = FALSE) {
  if (!size_col %in% names(plot_df) || !any(is.finite(suppressWarnings(as.numeric(plot_df[[size_col]]))))) {
    add_missing(stem, paste("no finite values for point size:", size_col))
    return(invisible(NULL))
  }
  d <- plot_df %>%
    filter(is.finite(log2_or_plot), is.finite(rank_score), !is.na(tail_plot), tail_plot %in% names(tail_shapes)) %>%
    mutate(size_value = suppressWarnings(as.numeric(.data[[size_col]]))) %>%
    filter(is.finite(size_value))
  if (abs_size) d <- d %>% mutate(size_value = abs(size_value))
  d <- d %>%
    arrange(namespace_plot, tail_plot, log2_or_plot) %>%
    mutate(
      term_tail_label = paste0(term_label_plain, " [", tail_plot, "]"),
      term_tail_label = factor(term_tail_label, levels = rev(unique(term_tail_label)))
    )

  bar_spec <- tai_bar_spec(size_col)
  has_bar <- !is.null(bar_spec) && any(is.finite(suppressWarnings(as.numeric(d[[bar_spec$column]]))))
  if (has_bar) {
    overlay <- prepare_tai_bar_overlay(d, bar_spec$column)
    d <- overlay$data
  }

  p <- ggplot(d, aes(y = term_tail_label))
  if (has_bar) {
    p <- p + geom_col(
      aes(x = bar_scaled), orientation = "y", width = 0.72,
      fill = "#BFDFF2", colour = "#6BAED6", linewidth = 0.18, alpha = 0.68
    )
  }
  p <- p +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.45) +
    geom_point(
      aes(x = log2_or_plot, size = size_value, shape = tail_plot, fill = rank_score),
      alpha = 0.94, colour = "grey15", stroke = 0.40
    ) +
    facet_grid(
      namespace_plot ~ ., scales = "free_y", space = "free_y", switch = "y",
      labeller = namespace_facet_labeller
    ) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, 38)) +
    scale_shape_manual(values = tail_shapes, drop = FALSE) +
    scale_fill_viridis_c(option = "inferno", direction = 1, name = math_labels("minus_log10_fdr")) +
    scale_size_continuous(range = size_range, name = size_title) +
    labs(
      y = NULL, shape = "tAI tail",
      title = "Selected GO terms from tAI-tail enrichment",
      subtitle = if (has_bar) {
        paste0(
          "Light-blue bars show ", tolower(bar_spec$title),
          " on the top axis; points show log2(CMH odds ratio), size and FDR support."
        )
      } else {
        "Shape identifies the tAI tail; fill shows FDR support; the dashed line denotes log2(odds ratio) = 0."
      }
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 13.0, lineheight = 1.02),
      axis.text.x = element_text(size = 13.5),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12.5),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 14.0, lineheight = 1.0),
      strip.background = element_rect(fill = "grey94", colour = "grey55", linewidth = 0.35),
      strip.placement = "outside",
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.box.just = "center",
      legend.title = element_text(size = 12.5, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.size = grid::unit(0.82, "cm"),
      plot.margin = margin(8, 12, 8, 12)
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(size = 5, fill = "grey75")),
      fill = guide_colourbar(order = 2, barwidth = grid::unit(5.0, "cm"), barheight = grid::unit(0.55, "cm")),
      size = guide_legend(order = 3, override.aes = list(shape = 21, fill = "grey75"))
    )
  if (has_bar) {
    p <- p + scale_x_continuous(
      name = math_labels("log2_or_cmh"),
      sec.axis = sec_axis(
        ~ . / overlay$scale_factor, name = bar_spec$title,
        breaks = overlay$breaks
      )
    )
  } else {
    p <- p + labs(x = math_labels("log2_or_cmh"))
  }
  save_plot_pair_local(p, stem, width = 17.5, height = max(8.5, 0.42 * nrow(d) + 3.8))
}

# Genome support is the least redundant size encoding and is used for the main
# figure. Additional requested versions use descriptive tail means.
if (any(is.finite(plot_df$support_plot))) {
  make_enrichment_dotplot("support_plot", "Informative genomes", "chosen_go_terms_enrichment_dotplot")
} else if ("n_tail_genomes" %in% names(plot_df)) {
  make_enrichment_dotplot("n_tail_genomes", "Tail genomes", "chosen_go_terms_enrichment_dotplot")
} else {
  make_enrichment_dotplot("rank_score", math_labels("minus_log10_fdr"), "chosen_go_terms_enrichment_dotplot")
}
if ("mean_tail_tAI" %in% names(plot_df)) {
  make_enrichment_dotplot("mean_tail_tAI", "Mean tAI", "chosen_go_terms_enrichment_dotplot_size_mean_tai", c(3.0, 10.0))
}
if ("mean_tail_tAI_z" %in% names(plot_df)) {
  make_enrichment_dotplot("mean_tail_tAI_z", "Absolute mean tAI z-score", "chosen_go_terms_enrichment_dotplot_size_mean_tai_z", c(3.0, 10.0), abs_size = TRUE)
}

# Overlay raw tAI bars and enrichment points in one panel. A secondary top axis
# preserves the raw-tAI units; bar lengths are linearly mapped only for display.
if (
  all(c("mean_tail_tAI", "mean_tail_tAI_z") %in% names(plot_df))
) {
  d <- plot_df %>%
    filter(
      is.finite(log2_or_plot), is.finite(rank_score),
      is.finite(mean_tail_tAI), is.finite(mean_tail_tAI_z),
      !is.na(tail_plot), tail_plot %in% names(tail_shapes)
    ) %>%
    arrange(namespace_plot, tail_plot, log2_or_plot) %>%
    mutate(
      term_tail_label = paste0(term_label_plain, " [", tail_plot, "]"),
      term_tail_label = factor(term_tail_label, levels = rev(unique(term_tail_label)))
    )

  if (nrow(d) > 0L) {
    overlay <- prepare_tai_bar_overlay(d, "mean_tail_tAI")
    d <- overlay$data
    combined_tai <- ggplot(d, aes(y = term_tail_label)) +
      geom_col(
        aes(x = bar_scaled), orientation = "y", width = 0.72,
        fill = "#BFDFF2", colour = "#6BAED6", linewidth = 0.18, alpha = 0.68
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.45) +
      geom_point(
        aes(
          x = log2_or_plot, size = abs(mean_tail_tAI_z),
          shape = tail_plot, fill = rank_score
        ),
        alpha = 0.94, colour = "grey15", stroke = 0.42
      ) +
      facet_grid(
        namespace_plot ~ ., scales = "free_y", space = "free_y", switch = "y",
        labeller = namespace_facet_labeller
      ) +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, 40)) +
      scale_x_continuous(
        name = math_labels("log2_or_cmh"),
        sec.axis = sec_axis(
          ~ . / overlay$scale_factor, name = "Mean raw tAI",
          breaks = overlay$breaks
        )
      ) +
      scale_shape_manual(values = tail_shapes, drop = FALSE) +
      scale_fill_viridis_c(option = "inferno", name = math_labels("minus_log10_fdr")) +
      scale_size_continuous(range = c(3.0, 10.0), name = "Absolute mean tAI z-score") +
      labs(
        y = NULL, shape = "tAI tail",
        title = "Selected GO terms: raw tAI and genome-stratified tail enrichment",
        subtitle = paste(
          "Light-blue bars show mean raw tAI on the top axis; points show log2(CMH odds ratio),",
          "with absolute mean tAI z-score encoded by size."
        )
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.y = element_text(size = 13.0, lineheight = 1.02),
        axis.text.x = element_text(size = 13.0),
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 14.0),
        strip.background = element_rect(fill = "grey94", colour = "grey55", linewidth = 0.35),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.title = element_text(size = 12.5, face = "bold"),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 12.5),
        plot.margin = margin(8, 12, 8, 12)
      )

    save_plot_pair_local(
      combined_tai,
      "chosen_go_terms_enrichment_raw_tai_bars_tai_z_dotplot",
      width = 18.5,
      height = max(9.5, 0.48 * nrow(d) + 4.4)
    )
  }
}

# A complementary encoding combines genome support and standardized tAI without
# adding another panel: point size is informative-genome count and fill is mean
# tail tAI_z. FDR remains a selection/ranking criterion in the input table.
if (
  all(c("mean_tail_tAI", "mean_tail_tAI_z", "support_plot") %in% names(plot_df)) &&
  any(is.finite(plot_df$mean_tail_tAI)) &&
  any(is.finite(plot_df$mean_tail_tAI_z)) && any(is.finite(plot_df$support_plot))
) {
  d <- plot_df %>%
    filter(
      is.finite(log2_or_plot), is.finite(mean_tail_tAI_z),
      is.finite(support_plot), !is.na(tail_plot), tail_plot %in% names(tail_shapes)
    ) %>%
    arrange(namespace_plot, tail_plot, log2_or_plot)
  d <- d %>% mutate(
    term_tail_label = paste0(term_label_plain, " [", tail_plot, "]"),
    term_tail_label = factor(term_tail_label, levels = rev(unique(term_tail_label)))
  )
  overlay <- prepare_tai_bar_overlay(d, "mean_tail_tAI")
  d <- overlay$data
  p_combined <- ggplot(d, aes(y = term_tail_label)) +
    geom_col(
      aes(x = bar_scaled), orientation = "y", width = 0.72,
      fill = "#BFDFF2", colour = "#6BAED6", linewidth = 0.18, alpha = 0.68
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.45) +
    geom_point(
      aes(x = log2_or_plot, size = support_plot, shape = tail_plot, fill = mean_tail_tAI_z),
      alpha = 0.94, colour = "grey15", stroke = 0.45
    ) +
    facet_grid(
      namespace_plot ~ ., scales = "free_y", space = "free_y", switch = "y",
      labeller = namespace_facet_labeller
    ) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, 38)) +
    scale_shape_manual(values = tail_shapes, drop = FALSE) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = "Mean tAI z-score"
    ) +
    scale_size_continuous(range = c(3.5, 10.5), name = "Informative genomes") +
    scale_x_continuous(
      name = math_labels("log2_or_cmh"),
      sec.axis = sec_axis(
        ~ . / overlay$scale_factor, name = "Mean raw tAI",
        breaks = overlay$breaks
      )
    ) +
    labs(
      y = NULL, shape = "tAI tail",
      title = "Selected GO terms: enrichment, genome support and tAI z-score",
      subtitle = "Light-blue bars show mean raw tAI on the top axis; point size shows informative genomes and fill shows mean within-genome standardized tAI."
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 13.0, lineheight = 1.02),
      axis.text.x = element_text(size = 13.5),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 14, lineheight = 1.0),
      strip.background = element_rect(fill = "grey94", colour = "grey55", linewidth = 0.35),
      strip.placement = "outside",
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.just = "center",
      legend.title = element_text(size = 12.5, face = "bold"),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12.5),
      plot.margin = margin(8, 12, 8, 12)
    )
  save_plot_pair_local(
    p_combined,
    "chosen_go_terms_enrichment_support_mean_tai_z",
    width = 9.0,
    height = max(14.0, 0.43 * nrow(d) + 3.8)
  )
}

if (any(is.finite(plot_df$support_plot))) {
  d <- plot_df %>% filter(is.finite(support_plot), !is.na(namespace_plot)) %>% arrange(namespace_plot, desc(support_plot))
  d$term_label_plain <- factor(d$term_label_plain, levels = rev(unique(d$term_label_plain)))
  p <- ggplot(d, aes(x = term_label_plain, y = support_plot, fill = namespace_plot)) +
    geom_col(colour = "grey30", linewidth = 0.15, width = 0.72) +
    coord_flip() +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, 46)) +
    scale_fill_manual(
      values = namespace_colours,
      labels = function(x) stringr::str_wrap(x, width = 12),
      na.value = "grey80"
    ) +
    labs(
      x = NULL, y = clean_label(support_col), fill = "GO namespace",
      title = "Genome support for selected GO terms"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 11.5, lineheight = 1.02),
      axis.text.x = element_text(size = 13),
      plot.title = element_text(size = 18, face = "bold"),
      legend.position = "bottom", legend.justification = "center", legend.box.just = "center",
      legend.title = element_text(size = 12.5, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.width = grid::unit(1.2, "cm"), plot.margin = margin(8, 12, 8, 12)
    ) +
    guides(fill = guide_legend(override.aes = list(colour = NA)))
  save_plot_pair_local(p, "chosen_go_terms_genome_support", width = 15.0, height = max(7.5, 0.38 * nrow(d) + 3.4))
} else {
  add_missing("chosen_go_terms_genome_support", "no finite genome-support column found")
}

make_gene_membership_plots <- function() {
  if (!file.exists(args$gene_table)) {
    add_missing("chosen_go_terms_tai_distributions", paste("gene table not found:", args$gene_table))
    return(invisible(NULL))
  }
  genes <- readr::read_tsv(args$gene_table, show_col_types = FALSE, progress = FALSE)
  required <- c("sample", "gene_id", "go_terms")
  missing_cols <- setdiff(required, names(genes))
  if (length(missing_cols) > 0L) {
    add_missing("chosen_go_terms_tai_distributions", paste("gene table lacks columns:", paste(missing_cols, collapse = ", ")))
    return(invisible(NULL))
  }

  value_cols <- intersect(c("tAI_z", "tAI"), names(genes))
  if (length(value_cols) == 0L) {
    add_missing("chosen_go_terms_tai_distributions", "gene table has neither tAI_z nor tAI")
    return(invisible(NULL))
  }

  # A GO term may be selected for both tails. Distribution plots represent all
  # genes annotated with the term, so each GO ID must appear only once here.
  chosen_small <- plot_df %>%
    distinct(go_id_plot, go_name_plot, namespace_plot, term_label_plain)

  go_long <- genes %>%
    transmute(
      sample = as.character(sample),
      gene_id = as.character(gene_id),
      go_id_plot = stringr::str_extract_all(toupper(as.character(go_terms)), "GO:[0-9]{7}"),
      tAI_z = if ("tAI_z" %in% names(genes)) suppressWarnings(as.numeric(.data[["tAI_z"]])) else NA_real_,
      tAI = if ("tAI" %in% names(genes)) suppressWarnings(as.numeric(.data[["tAI"]])) else NA_real_
    ) %>%
    tidyr::unnest_longer(go_id_plot) %>%
    filter(!is.na(go_id_plot), nzchar(go_id_plot)) %>%
    distinct(sample, gene_id, go_id_plot, .keep_all = TRUE) %>%
    inner_join(chosen_small, by = "go_id_plot")

  if (nrow(go_long) == 0L) {
    add_missing("chosen_go_terms_tai_distributions", "no genes matched selected GO terms")
    return(invisible(NULL))
  }

  counts <- go_long %>%
    group_by(go_id_plot, go_name_plot, namespace_plot, term_label_plain) %>%
    summarise(
      n_gene_term_memberships = n_distinct(paste(sample, gene_id, sep = "\t")),
      n_organisms = n_distinct(sample),
      .groups = "drop"
    ) %>%
    mutate(term_label_counts_plain = paste0(term_label_plain, " [genes=", n_gene_term_memberships, "; organisms=", n_organisms, "]"))

  readr::write_tsv(counts %>% mutate(across(where(is.character), clean_text)), file.path(args$output_dir, "chosen_go_terms_gene_membership_counts.tsv"), na = "NA")

  plot_long <- go_long %>%
    left_join(counts %>% select(go_id_plot, term_label_counts_plain, n_gene_term_memberships), by = "go_id_plot")

  make_distribution <- function(value_col, stem, y_lab) {
    df <- plot_long %>%
      filter(is.finite(.data[[value_col]])) %>%
      mutate(term_label_counts_plain = factor(term_label_counts_plain, levels = counts %>% arrange(namespace_plot, n_gene_term_memberships) %>% pull(term_label_counts_plain)))
    if (nrow(df) < 20L || n_distinct(df$term_label_counts_plain) < 2L) {
      add_missing(stem, paste("too few finite", value_col, "values"))
      return(invisible(NULL))
    }
    fill_values <- if (identical(value_col, "tAI")) tai_namespace_colours else namespace_colours
    p <- ggplot(df, aes(y = term_label_counts_plain, x = .data[[value_col]], fill = namespace_plot)) +
      geom_violin(
        orientation = "y", colour = "grey35", width = 0.72,
        trim = TRUE, scale = "width", alpha = 0.82
      ) +
      geom_boxplot(
        orientation = "y", width = 0.10, outlier.alpha = 0.15,
        fill = "white", colour = "grey25"
      ) +
      facet_grid(
        namespace_plot ~ ., scales = "free_y", space = "free_y", switch = "y",
        labeller = namespace_facet_labeller
      ) +
      scale_y_discrete(
        labels = function(x) stringr::str_wrap(x, 52),
        expand = expansion(add = 0.45)
      ) +
      scale_fill_manual(values = fill_values, na.value = "grey70", guide = "none") +
      labs(
        x = y_lab,
        y = NULL,
        title = if (value_col == "tAI_z") "Selected GO terms: tAI z-score distributions" else "Selected GO terms: tAI distributions"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        panel.spacing = grid::unit(0.7, "lines"),
        axis.text.y = element_text(size = 15.2, lineheight = 1.12),
        axis.text.x = element_text(size = 14.5),
        axis.title.x = element_text(size = 15.5),
        plot.title = element_text(size = 20, face = "bold"),
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 15.5, lineheight = 1.0),
        strip.background = element_rect(fill = "grey94", colour = "grey55", linewidth = 0.35),
        strip.placement = "outside",
        plot.margin = margin(10, 16, 10, 45)
      )
    save_plot_pair_local(p, stem, width = 22.0, height = max(12.0, 1.02 * n_distinct(df$term_label_counts_plain) + 5.5))
  }

  if ("tAI_z" %in% value_cols) {
    make_distribution("tAI_z", "chosen_go_terms_tai_z_distributions", metric_axis_label("tAI_z"))
  }
  if ("tAI" %in% value_cols) {
    make_distribution("tAI", "chosen_go_terms_tai_distributions", metric_axis_label("tAI"))
  }
}

make_gene_membership_plots()

if (nrow(missing) == 0L) missing <- tibble(plot = "none", reason = "all requested GO plots generated")
readr::write_tsv(missing, file.path(args$output_dir, "XXXXX_missing_go_plots.tsv"))

write_plot_manifest()
message("Done: ", args$output_dir)
