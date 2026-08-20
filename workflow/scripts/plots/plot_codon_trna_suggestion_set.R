#!/usr/bin/env Rscript
# Codon/tRNA exploratory plots reconstructed from legacy script suggestions.
# All measures and transformations are written to TSV diagnostics.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
})
source("workflow/scripts/lib/plot_style_helpers.R")

parse_args <- function() {
  x <- commandArgs(trailingOnly = TRUE)
  out <- list(codon_profiles = "results/tables/codon_profiles.tsv",
              genome_summary = "results/tables/genome_summary.tsv",
              output_dir = "results/plots/script_suggestions/codon_trna",
              formats = "png,pdf",
              trna_count_tables = character(),
              clean_trnascan_tables = character(),
              trna_audit_samples = character())
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Invalid CLI arguments")
    out[[gsub("-", "_", sub("^--", "", x[[i]]))]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}
args <- if (exists("snakemake")) {
  list(codon_profiles = snakemake@input[["codon_profiles"]],
       genome_summary = snakemake@input[["genome_summary"]],
       trna_count_tables = as.character(snakemake@input[["trna_count_tables"]] %||% character()),
       clean_trnascan_tables = as.character(snakemake@input[["clean_trnascan_tables"]] %||% character()),
       trna_audit_samples = as.character(snakemake@params[["trna_audit_samples"]] %||% character()),
       output_dir = snakemake@params[["output_dir"]] %||% "results/plots/script_suggestions/codon_trna",
       formats = paste(snakemake@params[["formats"]] %||% c("png", "pdf"), collapse = ","))
} else parse_args()

formats <- str_split(args$formats, "[,;[:space:]]+")[[1]] |> trimws()
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(args$output_dir, c("codon_variation_across_organisms_normalized.png",
                                    "codon_variation_across_organisms_normalized.pdf",
                                    "codon_variation_across_organisms_normalized.tsv",
                                    "codon_usage_per_species.png",
                                    "codon_usage_per_species.pdf",
                                    "codon_profiles_per_genome.png",
                                    "codon_profiles_per_genome.pdf")), force = TRUE)

save_multi <- function(p, stem, width, height) {
  p <- prepare_plot_for_export(p)
  for (fmt in formats[nzchar(formats)]) {
    ggsave(file.path(args$output_dir, paste0(stem, ".", fmt)), p,
           width = width, height = height, dpi = 300, limitsize = FALSE)
  }
}
pick_col <- function(d, candidates) {
  z <- intersect(candidates, names(d)); if (length(z)) z[[1]] else NA_character_
}
missing_plots <- tibble(plot = character(), reason = character())
skip <- function(plot, reason) {
  missing_plots <<- bind_rows(missing_plots, tibble(plot = plot, reason = reason))
  message("XXXXX skipped ", plot, ": ", reason)
}

if (!file.exists(args$codon_profiles)) stop("Missing codon_profiles table: ", args$codon_profiles)
codons <- read_tsv(args$codon_profiles, show_col_types = FALSE, progress = FALSE)
for (nm in c("sample", "codon")) if (!nm %in% names(codons)) stop("Missing column: ", nm)
for (nm in c("amino_acid", "phylum", "lifestyle", "trna_anticodon", "trna_id")) {
  if (!nm %in% names(codons)) codons[[nm]] <- NA_character_
}

genome_summary <- readr::read_tsv(args$genome_summary, show_col_types = FALSE, progress = FALSE)
if (!all(c("sample", "mean_tAI") %in% names(genome_summary))) {
  stop("genome_summary.tsv must contain sample and mean_tAI")
}
mean_tai_map <- genome_summary %>%
  transmute(sample = as.character(sample), mean_tAI = suppressWarnings(as.numeric(mean_tAI))) %>%
  distinct(sample, .keep_all = TRUE) %>%
  { stats::setNames(.$mean_tAI, .$sample) }

sample_lookup <- short_sample_labels(unique(codons$sample))
codons <- codons %>%
  mutate(
    codon = toupper(as.character(codon)),
    sample_label = unname(sample_lookup[as.character(sample)]),
    phylum = clean_label(phylum),
    lifestyle = clean_label(lifestyle),
    phylum = factor(phylum, levels = order_project_groups(phylum, "phylum"))
  )

# Use one display amino-acid family per codon. Alternative genetic codes are
# retained in a diagnostic table, while the modal annotation controls only the
# heatmap facet placement.
codon_amino_counts <- codons %>%
  filter(!is.na(codon), !is.na(amino_acid)) %>%
  count(codon, amino_acid, name = "n_rows") %>%
  arrange(codon, desc(n_rows), amino_acid)

ambiguous_map <- codon_amino_counts %>%
  group_by(codon) %>%
  mutate(n_amino_acids = n_distinct(amino_acid)) %>%
  filter(n_amino_acids > 1L) %>%
  ungroup()
if (nrow(ambiguous_map) > 0L) {
  readr::write_tsv(
    ambiguous_map,
    file.path(args$output_dir, "codon_amino_acid_assignment_diagnostics.tsv")
  )
}

codon_display_map <- codon_amino_counts %>%
  group_by(codon) %>%
  slice_head(n = 1L) %>%
  ungroup() %>%
  transmute(codon, amino_acid_plot = amino_acid)

codons <- codons %>%
  select(-amino_acid) %>%
  left_join(codon_display_map, by = "codon") %>%
  rename(amino_acid = amino_acid_plot)

trna_measure <- pick_col(codons, c("tRNA_weight", "trna_absolute_weight", "trna_copy_number"))
usage_measure <- pick_col(codons, c("codon_frequency", "codon_count", "genome_RSCU"))
readr::write_tsv(tibble(
  plot_family = c("tRNA heatmaps", "codon usage"),
  source_column = c(trna_measure, usage_measure),
  interpretation = c(
    "Codon-specific tRNA supply/weight. The standardized heatmap is a within-genome z-score across codons.",
    "codon_frequency = fraction of all CDS codons; codon_count = genome-wide count; genome_RSCU = relative synonymous codon usage."
  )), file.path(args$output_dir, "codon_trna_plot_methods.tsv"), na = "NA")

if (!is.na(trna_measure)) {
  trna_missing_diag <- codons %>%
    group_by(codon, amino_acid) %>%
    summarise(
      n_genomes_total = n_distinct(sample),
      n_genomes_finite = n_distinct(sample[is.finite(suppressWarnings(as.numeric(.data[[trna_measure]])))]),
      n_genomes_missing = n_genomes_total - n_genomes_finite,
      missing_fraction = n_genomes_missing / n_genomes_total,
      .groups = "drop"
    ) %>%
    arrange(desc(n_genomes_missing), codon)
  readr::write_tsv(
    trna_missing_diag,
    file.path(args$output_dir, "trna_measure_missingness_by_codon.tsv"),
    na = "NA"
  )
}

# ---- tRNA profiles ----------------------------------------------------------
if (!is.na(trna_measure)) {
  trna <- codons %>%
    mutate(value = suppressWarnings(as.numeric(.data[[trna_measure]]))) %>%
    filter(is.finite(value), !is.na(codon))
  if (nrow(trna)) {
    amino_levels <- sort(unique(as.character(trna$amino_acid)))
    trna <- trna %>% mutate(amino_acid = factor(amino_acid, levels = amino_levels))
    n_samples <- n_distinct(trna$sample)
    heat_height <- max(28.0, 6.5 + 0.26 * n_samples)

    heat_theme <- theme_minimal(base_size = 14) +
      theme(panel.grid = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14.0, face = "bold"),
            axis.text.y = element_text(size = 11.5),
            strip.text.x = element_text(face = "bold", size = 12.3),
            strip.text.y.left = element_text(angle = 90, face = "bold", size = 11.8, lineheight = 1.0),
            strip.background = element_rect(fill = "grey94", colour = "grey72"),
            strip.placement = "outside", legend.position = "bottom",
            panel.spacing.x = grid::unit(0.35, "lines"), panel.spacing.y = grid::unit(0.9, "lines"),
            legend.title = element_text(size = 12.5, face = "bold"),
            legend.text = element_text(size = 11.5),
            plot.title = element_text(size = 18, face = "bold"))

    message(
      "Raw tRNA/codon heatmap omitted here: the canonical codon_profiles/",
      "trna_weights_heatmap_by_phylum plot contains the same information."
    )

    trna_z <- trna %>% group_by(sample) %>%
      mutate(profile_sd = sd(value, na.rm = TRUE),
             value_z = if_else(is.finite(profile_sd) & profile_sd > 0,
                               (value - mean(value, na.rm = TRUE)) / profile_sd, NA_real_)) %>%
      ungroup()
    diag <- trna_z %>% distinct(sample, sample_label, profile_sd) %>%
      mutate(profile_constant = !is.finite(profile_sd) | profile_sd == 0)
    write_tsv(diag, file.path(args$output_dir, "codon_profiles_within_genome_standardization_diagnostics.tsv"), na = "NA")

    readr::write_tsv(
      tibble::tibble(
        output = "codon_profiles_per_genome_normalized",
        source_column = trna_measure,
        transformation = "Within each genome: (codon value - genome mean across codons) / genome SD across codons.",
        interpretation = "Positive values mark codons with above-genome-average tRNA supply; values are comparable as relative profile deviations, not raw copy numbers."
      ),
      file.path(args$output_dir, "codon_profiles_per_genome_normalized_method.tsv")
    )

    group_sizes <- trna_z %>%
      distinct(sample, phylum) %>%
      count(phylum, name = "n_genomes")
    keep_phyla <- group_sizes %>% filter(n_genomes >= 5L) %>% pull(phylum) %>% as.character()
    grouped <- trna_z %>%
      filter(is.finite(value_z)) %>%
      mutate(
        phylum_raw = as.character(phylum),
        phylum_group = if_else(phylum_raw %in% keep_phyla, phylum_raw, "Other"),
        mean_tAI = unname(mean_tai_map[as.character(sample)])
      )
    group_info <- grouped %>%
      distinct(sample, phylum_group) %>%
      count(phylum_group, name = "n_genomes") %>%
      mutate(
        phylum_group = factor(phylum_group, levels = order_project_groups(phylum_group, "phylum")),
        group_label = stringr::str_wrap(paste0(as.character(phylum_group), " (", n_genomes, ")"), 18)
      ) %>%
      arrange(phylum_group)
    group_label_map <- stats::setNames(group_info$group_label, as.character(group_info$phylum_group))
    sample_levels <- grouped %>%
      distinct(phylum_group, sample_label) %>%
      mutate(phylum_group = factor(phylum_group, levels = levels(group_info$phylum_group))) %>%
      arrange(phylum_group, sample_label) %>%
      pull(sample_label) %>%
      unique()
    codon_levels <- grouped %>%
      distinct(amino_acid, codon) %>%
      arrange(amino_acid, codon) %>%
      pull(codon) %>%
      unique()
    grouped <- grouped %>%
      mutate(
        phylum_group = factor(phylum_group, levels = levels(group_info$phylum_group)),
        group_label = factor(group_label_map[as.character(phylum_group)], levels = group_info$group_label),
        sample_label = factor(sample_label, levels = sample_levels),
        codon = factor(codon, levels = codon_levels)
      )
    readr::write_tsv(
      grouped %>% distinct(sample, sample_label, phylum_raw, phylum_group, group_label, mean_tAI),
      file.path(args$output_dir, "codon_profiles_per_genome_normalized_grouping.tsv"),
      na = "NA"
    )

    p_z_heat <- grouped %>%
      ggplot(aes(codon, sample_label, fill = value_z)) +
      geom_tile(width = 0.98, height = 0.98, colour = "white", linewidth = 0.035) +
      facet_grid(group_label ~ amino_acid, scales = "free", space = "free", switch = "y") +
      scale_fill_viridis_c(
        option = "C", begin = 0.12, end = 0.88, na.value = "#D9D9D9",
        guide = guide_colourbar(title.position = "top", barwidth = grid::unit(38, "cm"), barheight = grid::unit(0.65, "cm"))
      ) +
      scale_x_discrete(drop = TRUE, expand = expansion(add = 0)) + scale_y_discrete(drop = TRUE, expand = expansion(add = 0)) +
      labs(x = NULL, y = NULL, fill = "Within-genome z-score", title = NULL, subtitle = NULL) +
      heat_theme +
      theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y.left = element_blank(), strip.background.y = element_blank(),
        plot.margin = margin(8, 8, 8, 0)
      )

    p_z_tai <- grouped %>%
      distinct(sample, sample_label, group_label, mean_tAI) %>%
      ggplot(aes("Mean tAI", sample_label, fill = mean_tAI)) +
      geom_tile(width = 0.82, height = 0.98, colour = "white", linewidth = 0.08) +
      facet_grid(group_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
      scale_fill_viridis_c(
        option = "inferno", begin = 0.08, end = 0.92, na.value = "grey82",
        guide = guide_colourbar(title.position = "top", barwidth = grid::unit(7.5, "cm"), barheight = grid::unit(0.65, "cm"))
      ) +
      scale_x_discrete(position = "top", expand = expansion(add = 0.08)) +
      scale_y_discrete(drop = TRUE, expand = expansion(mult = 0)) +
      labs(x = NULL, y = NULL, fill = "Mean tAI") +
      heat_theme +
      theme(
        axis.text.x = element_text(angle = 0, size = 12.8, face = "bold"),
        strip.text.x = element_blank(), strip.background.x = element_blank(),
        plot.margin = margin(8, 4, 8, 8)
      )

    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Package 'patchwork' is required for the normalized codon heatmap.")
    }
    p_z <- (p_z_tai | p_z_heat) +
      patchwork::plot_layout(widths = c(4.8, 25), guides = "keep") +
      patchwork::plot_annotation(
        title = "Within-genome standardized tRNA/codon profile",
        subtitle = "Each cell is standardized within its genome across codons: (value - genome mean) / genome SD; phyla with fewer than five genomes are combined as Other."
      ) &
      theme(
        plot.title = element_text(size = 19, face = "bold"),
        plot.subtitle = element_text(size = 12.2),
        legend.position = "bottom"
      )
    save_multi(p_z, "codon_profiles_per_genome_normalized", 30.0, heat_height)

    variation <- trna %>% group_by(codon, amino_acid) %>%
      summarise(n_genomes = n_distinct(sample), mean_value = mean(value), sd_value = sd(value),
                cv = if_else(mean_value != 0, sd_value / abs(mean_value), NA_real_), .groups = "drop") %>%
      filter(is.finite(cv), n_genomes >= 3) %>% arrange(desc(cv))
    write_tsv(variation, file.path(args$output_dir, "codon_variation_across_organisms.tsv"), na = "NA")
    p_var <- variation %>% slice_head(n = 40L) %>%
      mutate(label = fct_reorder(paste0(codon, " (", amino_acid, ")"), cv)) %>%
      ggplot(aes(label, cv)) + geom_col(fill = "grey50", width = 0.72) + coord_flip() +
      labs(x = NULL, y = metric_cv_label(trna_measure),
           title = tex_label("Between-genome variability of codon-specific $\\mathrm{tRNA}$ profiles")) +
      theme_minimal(base_size = 11)
    save_multi(p_var, "codon_variation_across_organisms", 9.0, 8.5)
  } else skip("tRNA profile plots", paste("no finite", trna_measure, "values"))
} else skip("tRNA profile plots", "no tRNA weight/copy-number column")

# ---- Cys synonymous-codon usage on a common horizontal axis -----------------
if (!is.na(usage_measure)) {
  cys <- codons %>%
    filter(codon %in% c("TGT", "TGC")) %>%
    mutate(value = suppressWarnings(as.numeric(.data[[usage_measure]]))) %>%
    filter(is.finite(value))

  if (nrow(cys) >= 2 && n_distinct(cys$codon) == 2) {
    cys <- cys %>%
      group_by(sample) %>%
      mutate(
        cys_pair_total = sum(value, na.rm = TRUE),
        plot_value = if_else(cys_pair_total > 0, 100 * value / cys_pair_total, NA_real_)
      ) %>%
      ungroup() %>%
      filter(is.finite(plot_value))
    y_label <- "Share among Cys codons [%]"
    method_note <- paste0(
      "For each genome and codon c in {TGT,TGC}: 100 * ", usage_measure,
      "(c) / [", usage_measure, "(TGT) + ", usage_measure,
      "(TGC)]. For codon_count or codon_frequency this equals the percentage of Cys codons using c; for RSCU it is the synonymous-family share."
    )

    readr::write_tsv(
      tibble::tibble(plot = c(
                       "codon_usage_per_species_by_phylum",
                       "codon_usage_per_species_by_lifestyle",
                       "codon_usage_distribution"
                     ),
                      interpretation = method_note),
      file.path(args$output_dir, "cys_codon_usage_plot_interpretation.tsv")
    )
    readr::write_tsv(
      cys %>%
        transmute(
          sample, sample_label = as.character(sample_label), phylum = as.character(phylum),
          lifestyle, codon, source_measure = usage_measure, source_value = value,
          cys_pair_total, cys_codon_share_pct = plot_value
        ),
      file.path(args$output_dir, "cys_codon_usage_by_genome.tsv"),
      na = "NA"
    )

    make_species_plot <- function(group_col, stem, legend_title) {
      group_sizes <- cys %>%
        transmute(
          sample = as.character(sample),
          group_raw = dplyr::coalesce(clean_label(.data[[group_col]]), "Unknown")
        ) %>%
        distinct() %>%
        count(group_raw, name = "n_genomes")
      keep_groups <- group_sizes %>% filter(n_genomes >= 5L) %>% pull(group_raw)

      d <- cys %>%
        mutate(
          sample = as.character(sample),
          sample_label = as.character(sample_label),
          group_raw = dplyr::coalesce(clean_label(.data[[group_col]]), "Unknown"),
          group_value = if_else(group_raw %in% keep_groups, group_raw, "Other")
        )
      group_info <- d %>%
        distinct(sample, group_value) %>%
        count(group_value, name = "n_genomes") %>%
        mutate(
          group_value = factor(group_value, levels = order_project_groups(group_value, group_col)),
          group_label = stringr::str_wrap(
            paste0(as.character(group_value), " (", n_genomes, ")"),
            width = 18
          )
        ) %>%
        arrange(group_value)
      label_map <- stats::setNames(group_info$group_label, as.character(group_info$group_value))

      sample_order <- d %>%
        filter(codon == "TGC") %>%
        distinct(sample, sample_label, group_value, plot_value) %>%
        mutate(group_value = factor(group_value, levels = levels(group_info$group_value))) %>%
        arrange(group_value, plot_value, sample_label) %>%
        pull(sample_label) %>%
        unique()
      d <- d %>%
        mutate(
          group_value = factor(group_value, levels = levels(group_info$group_value)),
          group_label = factor(label_map[as.character(group_value)], levels = group_info$group_label),
          sample_label = factor(sample_label, levels = sample_order)
        )

      readr::write_tsv(
        d %>%
          distinct(sample, sample_label, original_group = group_raw, group = group_value, group_label) %>%
          arrange(group, sample_label),
        file.path(args$output_dir, paste0(stem, "_grouping_and_order.tsv"))
      )

      group_colours <- project_category_colours(d$group_value)
      p_species <- ggplot(
        d,
        aes(plot_value, sample_label, shape = codon, colour = group_value)
      ) +
        geom_line(aes(group = sample), linewidth = 0.42, alpha = 0.68) +
        geom_point(size = 2.8, alpha = 0.94) +
        facet_grid(group_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
        scale_shape_manual(values = c(TGT = 16, TGC = 17), drop = FALSE) +
        scale_colour_manual(values = group_colours, drop = FALSE) +
        scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
        labs(
          x = y_label, y = NULL, shape = "Cys codon", colour = legend_title,
          title = paste0("Cys synonymous-codon usage per genome, grouped by ", tolower(legend_title)),
          subtitle = paste0(
            "Codon usage is the within-genome percentage of Cys codons assigned to TGT or TGC, calculated from ",
            usage_measure,
            ". Groups with fewer than five genomes are combined as Other; genomes are ordered by TGC share."
          )
        ) +
        theme_minimal(base_size = 15) +
        theme(
          legend.position = "bottom",
          legend.box = "vertical",
          axis.text.x = element_text(size = 13.5),
          axis.text.y = element_text(size = 11.8),
          axis.title.x = element_text(size = 14.0),
          strip.text.y.left = element_text(angle = 90, face = "bold", size = 13.0, lineheight = 1.0),
          strip.background = element_rect(fill = "grey94", colour = "grey65", linewidth = 0.3),
          strip.placement = "outside",
          panel.spacing.y = grid::unit(0.9, "lines"),
          legend.title = element_text(size = 12.5, face = "bold"),
          legend.text = element_text(size = 12.0),
          plot.title = element_text(size = 17.5, face = "bold"),
          plot.subtitle = element_text(size = 12.2)
        )
      save_multi(p_species, stem, 16.0, max(9, 4.5 + 0.27 * n_distinct(d$sample)))
    }
    make_species_plot("phylum", "codon_usage_per_species_by_phylum", "Phylum")
    make_species_plot("lifestyle", "codon_usage_per_species_by_lifestyle", "Lifestyle")

    p_dist <- ggplot(cys, aes(codon, plot_value, fill = codon)) +
      geom_violin(trim = TRUE, alpha = 0.72) +
      geom_boxplot(width = 0.12, outlier.shape = NA) +
      scale_fill_manual(values = c(TGT = "#1F4E79", TGC = "#B8B8B8"), guide = "none") +
      scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
      labs(x = "Cys codon", y = y_label, title = "Distribution of Cys synonymous-codon usage across fungi", subtitle = NULL) +
      theme_minimal(base_size = 12)
    save_multi(p_dist, "codon_usage_distribution", 6.7, 5.4)
  } else skip("Cys codon usage", "TGT/TGC rows unavailable")
} else skip("Cys codon usage", "no codon-frequency/count/RSCU column")

# ---- Cys ACA/GCA anticodon copy number --------------------------------------
# ACA serves the Cys codon TGT and GCA serves TGC. Copy numbers are aggregated
# by anticodon, not by codon, so wobble mappings cannot be mistaken for genes.
make_horizontal_anticodon_plot <- function(data, value_col, point_col, title,
                                            x_label, stem, label_col = NULL,
                                            percent_axis = FALSE) {
  subtitle_text <- paste0(
    "Bars: ACA (serves TGT); points: GCA (serves TGC). ",
    "The log1p axis retains outlying genomes without a broken panel."
  )
  p <- ggplot(data, aes(x = .data[[value_col]], y = sample_label, fill = phylum)) +
    geom_col(width = 0.68, alpha = 0.90) +
    geom_point(
      aes(x = .data[[point_col]], shape = "GCA (serves TGC)"),
      size = 2.8, fill = "white", colour = "grey15", stroke = 0.7
    ) +
    scale_fill_manual(values = project_category_colours(data$phylum), drop = FALSE) +
    scale_shape_manual(values = c("GCA (serves TGC)" = 21), name = NULL) +
    scale_x_continuous(
      trans = "log1p",
      labels = if (percent_axis) function(x) paste0(signif(x, 3), "%") else scales::label_number(accuracy = 1),
      expand = expansion(mult = c(0, if (!is.null(label_col)) 0.28 else 0.15))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = x_label, y = NULL, fill = "Phylum", title = title,
      subtitle = subtitle_text
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom", legend.box = "vertical",
      axis.text.y = element_text(size = 8.4),
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 9.8, lineheight = 1.05),
      plot.margin = margin(7, if (!is.null(label_col)) 44 else 28, 7, 7)
    )
  if (!is.null(label_col)) {
    p <- p + geom_text(
      data = data,
      aes(x = .data[[value_col]], y = sample_label, label = .data[[label_col]]),
      inherit.aes = FALSE,
      hjust = -0.28,
      size = 2.7,
      colour = "grey20"
    )
  }
  save_multi(p, stem, 13.2, max(6.0, 2.8 + 0.30 * nrow(data)))
}

if ("trna_copy_number" %in% names(codons)) {
  copy_tbl <- codons %>%
    mutate(copy_number = suppressWarnings(as.numeric(trna_copy_number)))

  totals <- copy_tbl %>%
    filter(is.finite(copy_number), copy_number >= 0, !is.na(trna_id)) %>%
    group_by(sample, trna_id) %>%
    summarise(copy_number = max(copy_number, na.rm = TRUE), .groups = "drop") %>%
    group_by(sample) %>%
    summarise(total_trna_copies = sum(copy_number), .groups = "drop")

  cys_copy <- copy_tbl %>%
    mutate(trna_anticodon = toupper(as.character(trna_anticodon))) %>%
    filter(trna_anticodon %in% c("ACA", "GCA"), is.finite(copy_number), copy_number >= 0) %>%
    group_by(sample, sample_label, phylum, lifestyle, trna_anticodon) %>%
    summarise(anticodon_copy_number = max(copy_number, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = trna_anticodon,
      values_from = anticodon_copy_number,
      values_fill = 0
    ) %>%
    left_join(totals, by = "sample") %>%
    mutate(
      ACA = suppressWarnings(as.numeric(ACA)),
      GCA = suppressWarnings(as.numeric(GCA)),
      total_trna_copies = suppressWarnings(as.numeric(total_trna_copies)),
      cys_trna_copies = ACA + GCA,
      aca_percent_all_trna = 100 * ACA / total_trna_copies,
      gca_percent_all_trna = 100 * GCA / total_trna_copies,
      aca_percent_cys_trna = 100 * ACA / cys_trna_copies,
      copy_label = paste0(format(ACA, trim = TRUE), " / ", format(total_trna_copies, trim = TRUE))
    ) %>%
    filter(
      is.finite(ACA), is.finite(GCA), cys_trna_copies > 0,
      is.finite(total_trna_copies), total_trna_copies > 0
    )

  if (nrow(cys_copy)) {
    readr::write_tsv(
      cys_copy,
      file.path(args$output_dir, "cys_tgt_tgc_trna_copy_number_by_genome.tsv"),
      na = "NA"
    )

    cys_copy_with_tgt <- cys_copy %>% filter(ACA > 0)
    readr::write_tsv(
      cys_copy %>%
        filter(ACA <= 0) %>%
        select(sample, sample_label, phylum, lifestyle, ACA, GCA),
      file.path(args$output_dir, "cys_tgt_trna_excluded_genomes.tsv"),
      na = "NA"
    )

    if (nrow(cys_copy_with_tgt) > 0L) {
      count_plot <- cys_copy_with_tgt %>%
        arrange(ACA, sample_label) %>%
        mutate(sample_label = factor(sample_label, levels = unique(as.character(sample_label))))
      make_horizontal_anticodon_plot(
        count_plot,
        value_col = "ACA",
        point_col = "GCA",
        title = "Cys tRNA anticodon copy numbers per genome",
        x_label = "tRNA copy number",
        stem = "cys_tgt_trna_copy_number"
      )

      pct_plot <- cys_copy_with_tgt %>%
        arrange(aca_percent_all_trna, sample_label) %>%
        mutate(sample_label = factor(sample_label, levels = unique(as.character(sample_label))))
      make_horizontal_anticodon_plot(
        pct_plot,
        value_col = "aca_percent_all_trna",
        point_col = "gca_percent_all_trna",
        title = "Cys tRNA anticodons among elongator tRNAs used for tAI",
        x_label = "Anticodon copies among tAI-profile tRNA copies [%]",
        stem = "cys_tgt_trna_percent_of_all_trna",
        label_col = "copy_label",
        percent_axis = TRUE
      )

      composition <- cys_copy_with_tgt %>%
        select(sample, sample_label, phylum, lifestyle, ACA, GCA, aca_percent_cys_trna) %>%
        pivot_longer(c(ACA, GCA), names_to = "anticodon", values_to = "copy_number") %>%
        mutate(
          sample_label = factor(
            sample_label,
            levels = cys_copy_with_tgt %>% arrange(aca_percent_cys_trna, sample_label) %>% pull(sample_label) %>% as.character() %>% unique()
          ),
          anticodon = factor(anticodon, levels = c("GCA", "ACA"))
        )
      composition_plot <- ggplot(composition, aes(sample_label, copy_number, fill = anticodon)) +
        geom_col(position = "fill", width = 0.76) +
        coord_flip() +
        scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0))) +
        scale_fill_manual(
          values = c("GCA" = "#B8B8B8", "ACA" = "#1F4E79"),
          labels = c("GCA (serves TGC)", "ACA (serves TGT)")
        ) +
        labs(
          x = NULL, y = "Share within Cys tRNA copies", fill = "Cys anticodon",
          title = "Cys tRNA anticodon composition per genome"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.y = element_text(size = 8.2), legend.position = "bottom",
          plot.title = element_text(size = 15, face = "bold")
        )
      save_multi(composition_plot, "cys_trna_anticodon_composition", 11.5, max(6.0, 2.8 + 0.24 * nrow(cys_copy_with_tgt)))
    } else {
      skip("Cys TGT-serving tRNA plots", "no genome has a positive ACA anticodon copy number")
    }

    aca_fraction_median <- median(cys_copy$aca_percent_cys_trna, na.rm = TRUE)
    aca_fraction_plot <- ggplot(cys_copy, aes(x = aca_percent_cys_trna)) +
      geom_histogram(binwidth = 5, boundary = 0, fill = "#1F4E79", colour = "white", linewidth = 0.25) +
      geom_vline(xintercept = aca_fraction_median, linetype = "dashed", colour = "grey25", linewidth = 0.7) +
      scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
      labs(
        x = "ACA copies among Cys tRNA copies [%]", y = "Genomes",
        title = "Distribution of the tRNA-Cys(ACA) fraction",
        subtitle = paste0(
          "Denominator: ACA + GCA copies within each genome; dashed line: median = ",
          formatC(aca_fraction_median, digits = 3, format = "fg"), "%."
        )
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 13.0),
        plot.title = element_text(size = 17, face = "bold"),
        plot.subtitle = element_text(size = 11.8)
      )
    save_multi(aca_fraction_plot, "cys_trna_aca_fraction_distribution", 8.6, 6.2)

    readr::write_tsv(
      tibble::tibble(
        quantity = c("copy_label", "aca_percent_all_trna", "aca_percent_cys_trna"),
        definition = c(
          "ACA copy number / all elongator tRNA-gene copies represented in the tAI input profile; not ACA/GCA.",
          "100 * ACA / all elongator tRNA-gene copies represented in the tAI input profile.",
          "100 * ACA / (ACA + GCA); this is the within-Cys anticodon composition."
        )
      ),
      file.path(args$output_dir, "cys_trna_denominator_definitions.tsv")
    )
  } else skip("Cys anticodon copy plots", "no genome has a positive ACA or GCA copy number")
} else skip("Cys anticodon copy plots", "trna_copy_number column missing")

# ---- Focused tRNAscan-SE comparison for the two requested genomes -------------
# Counts come directly from cleaned tRNAscan-SE predictions. Cognate DNA codons
# are the reverse complements of anticodons; rows are ordered alphabetically by
# amino acid and then codon.
audit_samples <- as.character(args$trna_audit_samples %||% character())
clean_paths <- as.character(args$clean_trnascan_tables %||% character())

reverse_complement_codon <- function(anticodon) {
  anticodon <- as.character(anticodon)
  codon <- rep(NA_character_, length(anticodon))
  valid <- !is.na(anticodon) & grepl("^[ACGT]{3}$", anticodon)

  complement <- chartr("ACGT", "TGCA", anticodon[valid])
  codon[valid] <- paste0(
    substr(complement, 3L, 3L),
    substr(complement, 2L, 2L),
    substr(complement, 1L, 1L)
  )
  codon
}

format_plot_value <- function(x, digits) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

if (
  length(audit_samples) == 2L &&
  length(clean_paths) == length(audit_samples) &&
  all(file.exists(clean_paths))
) {
  raw_audit <- purrr::map2_dfr(clean_paths, audit_samples, function(path, sample_id) {
    tab <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
    needed <- c("trna_type", "anticodon", "score")
    if (!all(needed %in% names(tab))) {
      stop("Cleaned tRNAscan-SE table lacks trna_type/anticodon/score: ", path)
    }
    if (!"is_pseudo" %in% names(tab)) tab$is_pseudo <- FALSE

    tab %>%
      transmute(
        sample = sample_id,
        amino_acid = stringr::str_squish(as.character(trna_type)),
        anticodon = toupper(stringr::str_replace_all(
          stringr::str_squish(as.character(anticodon)), "U", "T"
        )),
        score = suppressWarnings(as.numeric(score)),
        is_pseudo = tolower(as.character(is_pseudo)) %in% c("true", "t", "1", "yes")
      ) %>%
      mutate(
        retained_prediction =
          !is_pseudo &
          !tolower(amino_acid) %in% c("imet", "fmet", "und") &
          stringr::str_detect(anticodon, "^[ACGT]{3}$"),
        codon = if_else(
          retained_prediction,
          reverse_complement_codon(anticodon),
          NA_character_
        )
      )
  })

  retained <- raw_audit %>%
    filter(retained_prediction, !is.na(codon)) %>%
    mutate(amino_acid = stringr::str_to_title(amino_acid))

  audit_summary <- raw_audit %>%
    group_by(sample) %>%
    summarise(
      n_cleaned_trnascan_predictions = n(),
      n_retained_trna_gene_predictions = sum(retained_prediction),
      mean_trnascan_score = if (
        any(retained_prediction & is.finite(score))
      ) mean(score[retained_prediction & is.finite(score)]) else NA_real_,
      .groups = "drop"
    ) %>%
    left_join(
      genome_summary %>%
        transmute(
          sample = as.character(sample),
          mean_tAI = suppressWarnings(as.numeric(mean_tAI))
        ),
      by = "sample"
    ) %>%
    mutate(
      organism = case_when(
        stringr::str_detect(sample, "^Zymoseptoria_tritici") ~ "Zymoseptoria tritici",
        stringr::str_detect(sample, "^Serpula_lacrymans") ~ "Serpula lacrymans",
        TRUE ~ clean_label(sample)
      ),
      facet_label = paste0(
        organism,
        "\nMean tAI = ", format_plot_value(mean_tAI, 3),
        "; mean tRNAscan-SE score = ", format_plot_value(mean_trnascan_score, 1)
      )
    )

  audit_subtitle <- audit_summary %>%
    arrange(match(sample, audit_samples)) %>%
    transmute(summary_label = paste0(
      organism,
      ": mean tRNAscan-SE score = ", format_plot_value(mean_trnascan_score, 1),
      "; mean tAI = ", format_plot_value(mean_tAI, 3)
    )) %>%
    pull(summary_label) %>%
    paste(collapse = "; ")

  catalogue <- retained %>%
    distinct(amino_acid, codon) %>%
    arrange(stringr::str_to_lower(amino_acid), codon) %>%
    mutate(
      codon_key = paste(amino_acid, codon, sep = "|"),
      codon_label = paste0(amino_acid, "\n", codon)
    )

  trna_counts <- retained %>%
    mutate(codon_key = paste(amino_acid, codon, sep = "|")) %>%
    count(sample, codon_key, name = "n_trna_genes")

  plot_data <- tidyr::crossing(
    sample = audit_samples,
    codon_key = catalogue$codon_key
  ) %>%
    left_join(catalogue, by = "codon_key") %>%
    left_join(trna_counts, by = c("sample", "codon_key")) %>%
    left_join(
      audit_summary %>% select(sample, organism, facet_label, mean_tAI, mean_trnascan_score),
      by = "sample"
    ) %>%
    mutate(
      n_trna_genes = dplyr::coalesce(n_trna_genes, 0L),
      codon_label = factor(codon_label, levels = catalogue$codon_label),
      organism = factor(
        organism,
        levels = audit_summary$organism[match(audit_samples, audit_summary$sample)]
      ),
      facet_label = factor(
        facet_label,
        levels = audit_summary$facet_label[match(audit_samples, audit_summary$sample)]
      )
    )

  readr::write_tsv(
    audit_summary,
    file.path(args$output_dir, "selected_genomes_trna_count_audit.tsv"),
    na = "NA"
  )
  readr::write_tsv(
    plot_data %>%
      transmute(
        sample, organism, amino_acid, codon,
        n_trna_genes, mean_tAI, mean_trnascan_score
      ),
    file.path(args$output_dir, "selected_genomes_trna_gene_copy_distribution.tsv"),
    na = "NA"
  )
  readr::write_tsv(
    tibble::tibble(
      field = c("n_trna_genes", "codon", "mean_tAI", "mean_trnascan_score"),
      interpretation = c(
        "Number of non-pseudogene elongator tRNA predictions in the cleaned tRNAscan-SE table.",
        "Cognate DNA codon obtained as the reverse complement of the predicted anticodon; wobble targets are not multiplied.",
        "Genome-wide mean gene tAI from genome_summary.tsv.",
        "Arithmetic mean tRNAscan-SE score across predictions counted in n_trna_genes."
      )
    ),
    file.path(args$output_dir, "selected_genomes_trna_count_audit_method.tsv")
  )
  unlink(
    file.path(args$output_dir, "selected_genomes_trna_prediction_class_counts.tsv"),
    force = TRUE
  )

  if (nrow(plot_data) > 0L) {
    organism_colours <- c(
      "Serpula lacrymans" = "#1F4E79",
      "Zymoseptoria tritici" = "#56B4E9"
    )
    p_audit <- ggplot(
      plot_data,
      aes(x = codon_label, y = n_trna_genes, fill = organism)
    ) +
      geom_col(
        position = position_dodge(width = 0.86),
        width = 0.78,
        colour = "grey25",
        linewidth = 0.12
      ) +
      scale_x_discrete(drop = FALSE, expand = expansion(add = 0.55)) +
      scale_y_continuous(
        breaks = scales::breaks_pretty(n = 6),
        expand = expansion(mult = c(0, 0.08))
      ) +
      scale_fill_manual(values = organism_colours, drop = FALSE) +
      labs(
        x = NULL,
        y = "tRNA genes predicted by tRNAscan-SE",
        fill = "Organism",
        title = "tRNA gene copy numbers in selected genomes",
        subtitle = paste0(
          audit_subtitle,
          ".\nCodons are reverse complements of predicted anticodons; initiator, undetermined and pseudogene records are excluded."
        )
      ) +
      theme_minimal(base_size = 17) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 17),
        legend.position = "bottom",
        legend.title = element_text(size = 15.5, face = "bold"),
        legend.text = element_text(size = 15),
        legend.key.width = grid::unit(0.9, "cm"),
        legend.spacing.x = grid::unit(0.35, "cm"),
        plot.title = element_text(size = 22, face = "bold"),
        plot.subtitle = element_text(size = 15, lineheight = 1.08)
      )
    save_multi(
      p_audit,
      "selected_genomes_trna_gene_copy_distribution",
      13.2,
      15.8
    )
  }
} else {
  skip(
    "selected_genomes_trna_gene_copy_distribution",
    "Serpula lacrymans and Zymoseptoria tritici cleaned tRNAscan-SE inputs were not both available"
  )
}

if (!nrow(missing_plots)) missing_plots <- tibble(plot = "none", reason = "all supported plots generated")
write_tsv(missing_plots, file.path(args$output_dir, "XXXXX_missing_plots.tsv"))
message("Done: ", args$output_dir)
