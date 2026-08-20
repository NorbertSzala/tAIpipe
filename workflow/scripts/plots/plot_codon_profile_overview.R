#!/usr/bin/env Rscript
# Codon/tRNA profile plots.
#
# Produces:
#   - codon_usage_variability.pdf/png
#   - large tRNA-weight heatmaps with codons on X and organisms on Y.
#
# Codon usage variability uses one explicitly selected usage column. By default
# the default is genome_RSCU. The selected column is written to
# codon_usage_variability_method.tsv so repeated runs are auditable.

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
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    codon_profiles = "results/tables/codon_profiles.tsv",
    genome_summary = "results/tables/genome_summary.tsv",
    output_dir = "results/plots/codon_profiles",
    formats = "png,pdf",
    heatmap_value = "auto",
    usage_value = "genome_RSCU",
    top_n_variable_codons = "30",
    heatmap_large_width = "16.5",
    heatmap_large_height = "23.4",
    heatmap_min_group_n = "5"
  )
  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--")) stop("Unexpected argument: ", args[[i]])
    key <- gsub("-", "_", sub("^--", "", args[[i]]))
    if (i == length(args)) stop("Missing value for argument: ", args[[i]])

    # Snakemake may render a Python list as either
    #   --formats png pdf
    # or, with quoting/formatting differences,
    #   --formats png,pdf / --formats pngpdf.
    # Accept all three forms so the plotting script is robust to the caller.
    if (identical(key, "formats")) {
      j <- i + 1L
      vals <- character()
      while (j <= length(args) && !startsWith(args[[j]], "--")) {
        vals <- c(vals, args[[j]])
        j <- j + 1L
      }
      if (length(vals) == 0L) stop("Missing value for argument: ", args[[i]])
      out[[key]] <- paste(vals, collapse = ",")
      i <- j
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

normalise_formats <- function(x) {
  x <- paste(as.character(x), collapse = ",")
  x <- trimws(x)
  if (identical(x, "pngpdf")) x <- "png,pdf"
  if (identical(x, "pdfpng")) x <- "pdf,png"
  out <- unlist(strsplit(x, "[,;[:space:]]+", perl = TRUE), use.names = FALSE)
  out <- trimws(out)
  out <- out[nzchar(out)]
  out <- unique(tolower(out))
  bad <- setdiff(out, c("png", "pdf", "svg"))
  if (length(bad) > 0L) stop("Unsupported output format(s): ", paste(bad, collapse = ", "))
  if (length(out) == 0L) out <- c("png", "pdf")
  out
}

args <- if (exists("snakemake")) {
  list(
    codon_profiles = snakemake@input[["codon_profiles"]] %||% snakemake@input[[1]],
    genome_summary = snakemake@input[["genome_summary"]] %||% "results/tables/genome_summary.tsv",
    output_dir = snakemake@params[["output_dir"]] %||% "results/plots/codon_profiles",
    formats = paste(snakemake@params[["formats"]] %||% c("png", "pdf"), collapse = ","),
    heatmap_value = snakemake@params[["heatmap_value"]] %||% "auto",
    usage_value = snakemake@params[["usage_value"]] %||% "genome_RSCU",
    top_n_variable_codons = as.character(snakemake@params[["top_n_variable_codons"]] %||% 30),
    heatmap_large_width = as.character(snakemake@params[["heatmap_large_width"]] %||% 16.5),
    heatmap_large_height = as.character(snakemake@params[["heatmap_large_height"]] %||% 23.4),
    heatmap_min_group_n = as.character(snakemake@params[["heatmap_min_group_n"]] %||% 5)
  )
} else parse_args()

formats <- normalise_formats(args$formats)
outdir <- args$output_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# Remove obsolete pre-revision clustered heatmaps so they are not mistaken for
# the tracked dendrogram outputs produced below.
unlink(list.files(
  outdir,
  pattern = "^trna_weights_heatmap_by_(phylum|lifestyle)_clustered",
  full.names = TRUE
), force = TRUE)

top_n_variable_codons <- as.integer(args$top_n_variable_codons)
heatmap_large_width <- as.numeric(args$heatmap_large_width)
heatmap_large_height <- as.numeric(args$heatmap_large_height)
heatmap_min_group_n <- as.integer(args$heatmap_min_group_n)

standard_code <- c(
  TTT="Phe", TTC="Phe", TTA="Leu", TTG="Leu",
  TCT="Ser", TCC="Ser", TCA="Ser", TCG="Ser",
  TAT="Tyr", TAC="Tyr", TAA="Stop", TAG="Stop",
  TGT="Cys", TGC="Cys", TGA="Stop", TGG="Trp",
  CTT="Leu", CTC="Leu", CTA="Leu", CTG="Leu",
  CCT="Pro", CCC="Pro", CCA="Pro", CCG="Pro",
  CAT="His", CAC="His", CAA="Gln", CAG="Gln",
  CGT="Arg", CGC="Arg", CGA="Arg", CGG="Arg",
  ATT="Ile", ATC="Ile", ATA="Ile", ATG="Met",
  ACT="Thr", ACC="Thr", ACA="Thr", ACG="Thr",
  AAT="Asn", AAC="Asn", AAA="Lys", AAG="Lys",
  AGT="Ser", AGC="Ser", AGA="Arg", AGG="Arg",
  GTT="Val", GTC="Val", GTA="Val", GTG="Val",
  GCT="Ala", GCC="Ala", GCA="Ala", GCG="Ala",
  GAT="Asp", GAC="Asp", GAA="Glu", GAG="Glu",
  GGT="Gly", GGC="Gly", GGA="Gly", GGG="Gly"
)

# Alphabetical amino-acid order; synonymous codons are alphabetical within amino acid.
amino_order <- sort(unique(unname(standard_code[standard_code != "Stop"])))

choose_col <- function(data, candidates, requested = "auto", label) {
  if (!identical(requested, "auto")) {
    if (!requested %in% names(data)) stop(label, " column not found: ", requested)
    return(requested)
  }
  hit <- candidates[candidates %in% names(data)][1]
  if (is.na(hit)) stop(label, " column not found. Tried: ", paste(candidates, collapse = ", "))
  hit
}

profiles <- readr::read_tsv(args$codon_profiles, show_col_types = FALSE, progress = FALSE)
required <- c("sample", "codon")
missing <- setdiff(required, names(profiles))
if (length(missing) > 0L) stop("codon_profiles.tsv lacks columns: ", paste(missing, collapse = ", "))

genome_summary <- readr::read_tsv(args$genome_summary, show_col_types = FALSE, progress = FALSE)
missing_genome <- setdiff(c("sample", "mean_tAI"), names(genome_summary))
if (length(missing_genome) > 0L) {
  stop("genome_summary.tsv lacks columns needed for the row annotation: ", paste(missing_genome, collapse = ", "))
}
mean_tai_lookup <- genome_summary %>%
  transmute(sample = as.character(sample), mean_tAI = suppressWarnings(as.numeric(mean_tAI))) %>%
  distinct(sample, .keep_all = TRUE)
mean_tai_map <- stats::setNames(mean_tai_lookup$mean_tAI, mean_tai_lookup$sample)

profiles <- profiles %>%
  mutate(
    codon = toupper(as.character(codon)),
    amino_acid = if ("amino_acid" %in% names(.)) as.character(amino_acid) else unname(standard_code[codon]),
    amino_acid = if_else(is.na(amino_acid), "Unknown", amino_acid),
    phylum = if ("phylum" %in% names(.)) clean_label(phylum) else "Unknown phylum",
    lifestyle = if ("lifestyle" %in% names(.)) clean_label(lifestyle) else "Unknown lifestyle"
  ) %>%
  filter(amino_acid != "Stop")

# A single codon may be annotated with different amino acids in genomes using
# alternative genetic codes. For a compact heatmap every codon must occupy one
# display facet, so use its most frequent annotation and write all conflicts to
# a diagnostic table. This changes display grouping only, never the values.
codon_amino_counts <- profiles %>%
  count(codon, amino_acid, name = "n_rows") %>%
  arrange(codon, desc(n_rows), amino_acid)

ambiguous_codon_map <- codon_amino_counts %>%
  group_by(codon) %>%
  mutate(n_amino_acids = n_distinct(amino_acid)) %>%
  filter(n_amino_acids > 1L) %>%
  ungroup()

if (nrow(ambiguous_codon_map) > 0L) {
  readr::write_tsv(
    ambiguous_codon_map,
    file.path(outdir, "trna_weights_ambiguous_codon_amino_acid_map.tsv")
  )
}

codon_display_map <- codon_amino_counts %>%
  group_by(codon) %>%
  slice_head(n = 1L) %>%
  ungroup() %>%
  transmute(codon, amino_acid_plot = amino_acid)

profiles <- profiles %>%
  select(-amino_acid) %>%
  left_join(codon_display_map, by = "codon") %>%
  rename(amino_acid = amino_acid_plot)

amino_order <- c(
  amino_order,
  sort(setdiff(unique(profiles$amino_acid), amino_order))
)

heatmap_col <- choose_col(
  profiles,
  c("tRNA_weight", "trna_weight", "trna_absolute_weight", "trna_copy_number", "relative_trna"),
  args$heatmap_value,
  "Heatmap value"
)
usage_col <- choose_col(
  profiles,
  c("genome_RSCU", "RSCU", "rscu", "codon_frequency", "frequency", "codon_count"),
  args$usage_value,
  "Codon usage variability value"
)

profiles <- profiles %>%
  mutate(
    heatmap_value = suppressWarnings(as.numeric(.data[[heatmap_col]])),
    usage_value = suppressWarnings(as.numeric(.data[[usage_col]]))
  )

readr::write_tsv(
  tibble(
    heatmap_value_column = heatmap_col,
    codon_usage_variability_column = usage_col,
    interpretation = ifelse(grepl("RSCU|rscu", usage_col),
      "CV of RSCU across genomes: variability of synonymous codon preference.",
      "CV of codon frequency/count across genomes: variability of raw codon abundance, more composition-dependent."
    ),
    heatmap_interpretation = "tRNA weights are supply-side codon-adaptation parameters used by tAI; they are not gene-level tAI and not codon usage like RSCU."
  ),
  file.path(outdir, "codon_usage_variability_method.tsv")
)

# Main codon-usage variability figure.
cv_tbl <- profiles %>%
  filter(is.finite(usage_value)) %>%
  group_by(codon, amino_acid) %>%
  summarise(
    mean_usage = mean(usage_value, na.rm = TRUE),
    sd_usage = sd(usage_value, na.rm = TRUE),
    cv = sd_usage / abs(mean_usage),
    n_genomes = n_distinct(sample),
    .groups = "drop"
  ) %>%
  filter(is.finite(cv), n_genomes >= 3L) %>%
  arrange(desc(cv)) %>%
  mutate(codon_label = paste0(codon, " (", amino_acid, ")"))

readr::write_tsv(cv_tbl, file.path(outdir, "codon_usage_variability.tsv"))

p_cv <- cv_tbl %>%
  slice_head(n = min(top_n_variable_codons, nrow(cv_tbl))) %>%
  mutate(codon_label = forcats::fct_reorder(codon_label, cv)) %>%
  ggplot(aes(x = codon_label, y = cv)) +
  geom_col(width = 0.75, fill = "#666666") +
  coord_flip() +
  labs(
    x = NULL,
    y = if (identical(usage_col, "genome_RSCU")) {
      "Codon variation based on genomic RSCU"
    } else {
      paste("Codon variation based on", metric_plain_label(usage_col))
    },
    title = "Most variable codons across genomes",
    subtitle = NULL
  ) +
  theme_minimal(base_size = 11)

save_plot_pair(p_cv, "codon_usage_variability", outdir, 7.8, 6.4, formats)

prepare_heatmap_df <- function(group_col, min_count = NULL) {
  group_name <- rlang::as_string(rlang::ensym(group_col))
  label_lookup <- short_sample_labels(unique(profiles$sample))
  df <- profiles %>%
    filter(is.finite(heatmap_value)) %>%
    mutate(
      group_raw = clean_label(.data[[group_name]]),
      sample_short = unname(label_lookup[as.character(sample)]),
      mean_tAI = unname(mean_tai_map[as.character(sample)]),
      amino_acid = factor(amino_acid, levels = amino_order),
      codon_label = codon
    )

  other_members <- character()
  if (!is.null(min_count)) {
    tab <- df %>% distinct(sample, group_raw) %>% count(group_raw, name = "n")
    keep <- tab %>% filter(n >= min_count) %>% pull(group_raw)
    other_members <- tab %>% filter(!group_raw %in% keep) %>% pull(group_raw)
    df <- df %>% mutate(group_raw = if_else(group_raw %in% keep, group_raw, "Other"))
  }

  group_labels <- df %>%
    distinct(sample, group_raw) %>%
    count(group_raw, name = "n") %>%
    mutate(
      group_label = stringr::str_wrap(paste0(group_raw, " (", n, ")"), width = 18),
      group_order = match(group_raw, order_project_groups(group_raw, group_name))
    ) %>%
    arrange(group_order, group_label)
  label_map <- stats::setNames(group_labels$group_label, group_labels$group_raw)

  codon_levels <- df %>% distinct(amino_acid, codon) %>% arrange(amino_acid, codon) %>% pull(codon)
  sample_levels <- df %>%
    distinct(group_raw, sample_short) %>%
    mutate(group_raw = factor(group_raw, levels = group_labels$group_raw)) %>%
    arrange(group_raw, sample_short) %>%
    pull(sample_short)

  out <- df %>%
    mutate(
      group_label = factor(label_map[group_raw], levels = group_labels$group_label),
      sample_short = factor(sample_short, levels = unique(sample_levels)),
      codon_label = factor(codon_label, levels = unique(codon_levels)),
      amino_acid = factor(amino_acid, levels = amino_order)
    )
  attr(out, "other_members") <- other_members
  out
}

heatmap_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "grey72", fill = NA, linewidth = 0.25),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 24, face = "bold"),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_blank(),
      strip.text.x = element_text(size = 18, face = "bold"),
      strip.text.y.left = element_text(size = 16, face = "bold", angle = 90, hjust = 0.5, lineheight = 1.0),
      strip.background = element_rect(fill = "grey94", colour = "grey70", linewidth = 0.25),
      strip.placement = "outside",
      legend.position = "bottom",
      panel.spacing.x = grid::unit(0.38, "lines"),
      panel.spacing.y = grid::unit(0.85, "lines"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(16),
      plot.title = element_text(size = 26, face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

make_large_heatmap <- function(group_col, suffix, min_count = NULL) {
  group_name <- rlang::as_string(rlang::ensym(group_col))
  df <- prepare_heatmap_df({{ group_col }}, min_count = min_count)
  other <- attr(df, "other_members")
  if (length(other) > 0L) {
    readr::write_tsv(
      tibble::tibble(group_variable = group_name, collapsed_to = "Other", original_group = clean_label(other)),
      file.path(outdir, paste0("other_members_by_", group_name, suffix, ".tsv"))
    )
  }

  readr::write_tsv(
    df %>% distinct(sample, sample_short, group_label) %>% arrange(group_label, sample_short),
    file.path(outdir, paste0("sample_order_by_", group_name, suffix, ".tsv"))
  )

  p_heat <- ggplot(df, aes(x = codon_label, y = sample_short, fill = heatmap_value)) +
    geom_tile(width = 0.98, height = 0.98, colour = "white", linewidth = 0.035) +
    facet_grid(group_label ~ amino_acid, scales = "free", space = "free", switch = "y") +
    scale_fill_viridis_c(
      option = "C", begin = 0.12, end = 0.88, na.value = "#D9D9D9",
      guide = guide_colourbar(
        order = 1,
        title.position = "top",
        barwidth = grid::unit(38, "cm"),
        barheight = grid::unit(0.65, "cm")
      )
    ) +
    scale_x_discrete(drop = TRUE, expand = expansion(add = 0)) +
    scale_y_discrete(drop = TRUE, expand = expansion(mult = 0)) +
    labs(
      x = NULL, y = NULL, fill = metric_axis_label(heatmap_col),
      title = NULL, subtitle = NULL
    ) +
    heatmap_theme() +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      strip.text.y.left = element_blank(), strip.background.y = element_blank(),
      plot.margin = margin(8, 8, 8, 0)
    )

  p_tai <- df %>%
    distinct(sample, sample_short, group_label, mean_tAI) %>%
    ggplot(aes(x = "Mean tAI", y = sample_short, fill = mean_tAI)) +
    geom_tile(width = 0.82, height = 0.98, colour = "white", linewidth = 0.08) +
    facet_grid(group_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_viridis_c(
      option = "D", begin = 0.08, end = 0.92, na.value = "grey82",
      guide = guide_colourbar(
        order = 2,
        title.position = "top",
        barwidth = grid::unit(38, "cm"),
        barheight = grid::unit(0.65, "cm")
      )
    ) +
    scale_x_discrete(position = "top", expand = expansion(add = 0.08)) +
    scale_y_discrete(drop = TRUE, expand = expansion(mult = 0)) +
    labs(x = NULL, y = NULL, fill = "Mean tAI") +
    heatmap_theme() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold"),
      axis.text.y = element_text(size = 18),
      strip.text.x = element_blank(), strip.background.x = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(8, 4, 8, 8)
    )

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for the tAI-annotated heatmaps.")
  }
  p <- (p_tai | p_heat) +
    patchwork::plot_layout(widths = c(6.5, 25), guides = "collect") +
    patchwork::plot_annotation() &
    theme(
      plot.title = element_text(size = 28, face = "bold"),
      legend.position = "bottom", legend.box = "vertical"
    )

  stem <- paste0("trna_weights_heatmap_by_", group_name, suffix, "_large_codon_x")
  # Keep individual cells legible on large fungal panels. The dynamic height
  # grows with the number of genomes instead of compressing 183 rows into a
  # fixed-height panel.
  dynamic_height <- max(heatmap_large_height, 6.5 + 0.24 * dplyr::n_distinct(df$sample_short))
  dynamic_width <- max(heatmap_large_width, 30.0)
  save_plot_pair(p, stem, outdir, dynamic_width, dynamic_height, formats)
}

make_large_heatmap(phylum, "")
make_large_heatmap(lifestyle, "")
make_large_heatmap(phylum, paste0("_count", heatmap_min_group_n), min_count = heatmap_min_group_n)
make_large_heatmap(lifestyle, paste0("_count", heatmap_min_group_n), min_count = heatmap_min_group_n)

hclust_segments <- function(hc) {
  n <- length(hc$order)
  leaf_y <- numeric(n)
  leaf_y[hc$order] <- seq_len(n)
  node_y <- numeric(n - 1L)
  rows <- vector("list", 3L * (n - 1L))
  k <- 1L
  child_info <- function(id) {
    if (id < 0L) list(y = leaf_y[-id], height = 0) else list(y = node_y[id], height = hc$height[id])
  }
  for (i in seq_len(n - 1L)) {
    left <- child_info(hc$merge[i, 1])
    right <- child_info(hc$merge[i, 2])
    parent_h <- hc$height[i]
    node_y[i] <- mean(c(left$y, right$y))
    rows[[k]] <- tibble(x = left$height, xend = parent_h, y = left$y, yend = left$y); k <- k + 1L
    rows[[k]] <- tibble(x = right$height, xend = parent_h, y = right$y, yend = right$y); k <- k + 1L
    rows[[k]] <- tibble(x = parent_h, xend = parent_h, y = left$y, yend = right$y); k <- k + 1L
  }
  bind_rows(rows)
}

cluster_association_results <- tibble()

make_dendrogram_heatmap <- function(annotation_col) {
  annotation_name <- rlang::as_string(rlang::ensym(annotation_col))
  sample_lookup <- short_sample_labels(unique(profiles$sample))
  base <- profiles %>%
    filter(is.finite(heatmap_value)) %>%
    mutate(
      sample_short = unname(sample_lookup[as.character(sample)]),
      annotation = clean_label(.data[[annotation_name]])
    ) %>%
    group_by(sample, sample_short, annotation, codon, amino_acid) %>%
    summarise(heatmap_value = mean(heatmap_value, na.rm = TRUE), .groups = "drop")

  wide <- base %>% select(sample, codon, heatmap_value) %>% pivot_wider(names_from = codon, values_from = heatmap_value)
  mat <- as.data.frame(wide)
  rn <- mat$sample
  mat$sample <- NULL
  mat <- as.matrix(mat)
  rownames(mat) <- rn
  for (j in seq_len(ncol(mat))) {
    med <- median(mat[, j], na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    mat[!is.finite(mat[, j]), j] <- med
  }
  if (nrow(mat) < 3L) return(invisible(NULL))

  # Spearman correlation is undefined for a genome with a constant codon
  # profile. Remove such rows before clustering and report them explicitly.
  row_variable <- apply(mat, 1L, function(x) {
    x <- x[is.finite(x)]
    length(x) >= 2L && length(unique(x)) >= 2L
  })

  omitted_samples <- rownames(mat)[!row_variable]
  if (length(omitted_samples) > 0L) {
    readr::write_tsv(
      tibble::tibble(
        sample = omitted_samples,
        reason = "constant_or_noninformative_codon_weight_profile"
      ),
      file.path(
        outdir,
        paste0(
          "trna_weights_dendrogram_omitted_samples_by_",
          annotation_name,
          ".tsv"
        )
      )
    )
  }

  mat <- mat[row_variable, , drop = FALSE]
  if (nrow(mat) < 3L) {
    stop(
      "Fewer than three informative genome profiles remain for ",
      annotation_name,
      " after removing constant profiles."
    )
  }

  profile_cor_raw <- suppressWarnings(
    stats::cor(
      t(mat),
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )

  n_profiles <- nrow(mat)
  if (length(profile_cor_raw) != n_profiles * n_profiles) {
    stop(
      "Unexpected correlation result for ", annotation_name,
      ": mat=", paste(dim(mat), collapse = "x"),
      "; correlation length=", length(profile_cor_raw),
      "; expected=", n_profiles * n_profiles
    )
  }

  # Reconstruct dimensions explicitly. This prevents accidental loss of matrix
  # attributes before conversion to a dist object.
  profile_cor <- matrix(
    as.numeric(profile_cor_raw),
    nrow = n_profiles,
    ncol = n_profiles,
    dimnames = list(rownames(mat), rownames(mat))
  )
  profile_cor[!is.finite(profile_cor)] <- 0
  profile_cor <- (profile_cor + t(profile_cor)) / 2
  diag(profile_cor) <- 1

  profile_dist <- matrix(
    1 - as.numeric(profile_cor),
    nrow = n_profiles,
    ncol = n_profiles,
    dimnames = dimnames(profile_cor)
  )
  profile_dist[!is.finite(profile_dist)] <- 1
  profile_dist[profile_dist < 0] <- 0
  profile_dist[profile_dist > 2] <- 2
  profile_dist <- (profile_dist + t(profile_dist)) / 2
  diag(profile_dist) <- 0

  if (
    !is.matrix(profile_dist) ||
    nrow(profile_dist) != ncol(profile_dist) ||
    nrow(profile_dist) != n_profiles
  ) {
    stop(
      "Internal error: non-square profile distance matrix for ",
      annotation_name,
      ". Dimensions: ",
      paste(dim(profile_dist), collapse = "x")
    )
  }

  message(
    "Clustering ", annotation_name,
    ": mat=", paste(dim(mat), collapse = "x"),
    "; profile_cor=", paste(dim(profile_cor), collapse = "x"),
    "; profile_dist=", paste(dim(profile_dist), collapse = "x")
  )

  hc <- stats::hclust(
    stats::as.dist(profile_dist),
    method = "average"
  )

  # Quantify whether profiles are more similar within the same annotation group.
  # This is a label-permutation test on profile distances, not a phylogenetic test.
  meta_for_test <- base %>% distinct(sample, annotation) %>%
    filter(sample %in% rownames(profile_dist), !is.na(annotation), nzchar(annotation))
  idx <- match(meta_for_test$sample, rownames(profile_dist))
  d_sub <- profile_dist[idx, idx, drop = FALSE]
  labels_obs <- as.character(meta_for_test$annotation)
  upper <- upper.tri(d_sub)
  same_obs <- outer(labels_obs, labels_obs, `==`) & upper
  diff_obs <- outer(labels_obs, labels_obs, `!=`) & upper
  within_mean <- if (any(same_obs)) mean(d_sub[same_obs], na.rm = TRUE) else NA_real_
  between_mean <- if (any(diff_obs)) mean(d_sub[diff_obs], na.rm = TRUE) else NA_real_
  separation <- between_mean - within_mean
  permutation_p <- NA_real_
  if (is.finite(separation) && length(unique(labels_obs)) >= 2L) {
    set.seed(1)
    perm_stats <- replicate(999L, {
      lab <- sample(labels_obs, replace = FALSE)
      same <- outer(lab, lab, `==`) & upper
      diff <- outer(lab, lab, `!=`) & upper
      mean(d_sub[diff], na.rm = TRUE) - mean(d_sub[same], na.rm = TRUE)
    })
    permutation_p <- (1 + sum(perm_stats >= separation, na.rm = TRUE)) / (1 + sum(is.finite(perm_stats)))
  }
  cluster_association_results <<- bind_rows(
    cluster_association_results,
    tibble(
      annotation = annotation_name,
      n_genomes = length(labels_obs),
      n_groups = n_distinct(labels_obs),
      mean_within_group_distance = within_mean,
      mean_between_group_distance = between_mean,
      between_minus_within = separation,
      between_to_within_ratio = between_mean / within_mean,
      permutation_p_value = permutation_p,
      permutations = 999L,
      interpretation = "Positive between-minus-within means profiles are more similar within annotation groups; this is not phylogenetic-tree concordance."
    )
  )

  ordered_samples <- rownames(mat)[hc$order]
  pos_map <- stats::setNames(seq_along(ordered_samples), ordered_samples)

  meta <- base %>%
    distinct(sample, sample_short, annotation) %>%
    filter(sample %in% ordered_samples) %>%
    mutate(
      y_pos = unname(pos_map[as.character(sample)]),
      mean_tAI = unname(mean_tai_map[as.character(sample)])
    ) %>%
    filter(is.finite(y_pos)) %>%
    arrange(y_pos)

  heat <- base %>%
    filter(sample %in% ordered_samples) %>%
    mutate(y_pos = unname(pos_map[as.character(sample)])) %>%
    filter(is.finite(y_pos))

  # A codon can have more than one amino-acid annotation across genomes,
  # for example because alternative genetic codes are represented. A single
  # heatmap x-axis nevertheless requires one display family per codon. Use the
  # most frequent annotation in the current dataset and record ambiguities.
  codon_amino_counts <- heat %>%
    count(codon, amino_acid, name = "n_rows") %>%
    arrange(codon, desc(n_rows), amino_acid)

  ambiguous_codon_map <- codon_amino_counts %>%
    group_by(codon) %>%
    mutate(n_amino_acids = n_distinct(amino_acid)) %>%
    filter(n_amino_acids > 1L) %>%
    ungroup()

  if (nrow(ambiguous_codon_map) > 0L) {
    readr::write_tsv(
      ambiguous_codon_map,
      file.path(
        outdir,
        paste0(
          "trna_weights_ambiguous_codon_amino_acid_map_by_",
          annotation_name,
          ".tsv"
        )
      )
    )
  }

  codon_display_map <- codon_amino_counts %>%
    group_by(codon) %>%
    slice_head(n = 1L) %>%
    ungroup() %>%
    transmute(codon, amino_acid_plot = amino_acid)

  amino_levels_plot <- c(
    amino_order,
    sort(setdiff(unique(codon_display_map$amino_acid_plot), amino_order))
  )

  heat <- heat %>%
    select(-amino_acid) %>%
    left_join(codon_display_map, by = "codon") %>%
    mutate(
      amino_acid = factor(amino_acid_plot, levels = amino_levels_plot)
    )

  codon_order <- heat %>%
    distinct(amino_acid, codon) %>%
    arrange(amino_acid, codon) %>%
    pull(codon) %>%
    unique()

  heat <- heat %>%
    mutate(codon = factor(codon, levels = codon_order))

  dend <- ggplot(hclust_segments(hc), aes(x = x, xend = xend, y = y, yend = yend)) +
    geom_segment(linewidth = 0.35, colour = "grey25") +
    scale_x_reverse(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0.5, nrow(meta) + 0.5), expand = expansion(mult = 0)) +
    labs(x = "1 - Spearman rho", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(5, 0, 5, 5))

  ann_colours <- project_category_colours(meta$annotation)
  ann <- ggplot(meta, aes(x = 1, y = y_pos, fill = annotation)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_manual(
      values = ann_colours,
      labels = function(x) stringr::str_wrap(x, width = 18),
      drop = FALSE,
      guide = guide_legend(order = 1, nrow = 3, byrow = TRUE)
    ) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0.5, nrow(meta) + 0.5), expand = expansion(mult = 0)) +
    labs(x = NULL, y = NULL, fill = clean_label(annotation_name)) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      plot.margin = margin(5, 0, 5, 0)
    )

  tai_ann <- ggplot(meta, aes(x = 1, y = y_pos, fill = mean_tAI)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_viridis_c(
      option = "D", begin = 0.08, end = 0.92, na.value = "grey82",
      guide = guide_colourbar(
        order = 3,
        title.position = "top",
        barwidth = grid::unit(38, "cm"),
        barheight = grid::unit(0.65, "cm")
      )
    ) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0.5, nrow(meta) + 0.5), expand = expansion(mult = 0)) +
    labs(x = NULL, y = NULL, fill = "Mean tAI", title = "Mean\ntAI") +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5, lineheight = 0.9),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      plot.margin = margin(5, 0, 5, 0)
    )

  label_map <- stats::setNames(meta$sample_short, meta$y_pos)
  hp <- ggplot(heat, aes(x = codon, y = y_pos, fill = heatmap_value)) +
    geom_tile(width = 0.98, height = 0.98, colour = "white", linewidth = 0.035) +
    facet_grid(. ~ amino_acid, scales = "free_x", space = "free_x") +
    scale_fill_viridis_c(
      option = "C", begin = 0.12, end = 0.88, na.value = "#D9D9D9",
      guide = guide_colourbar(
        order = 2,
        title.position = "top",
        barwidth = grid::unit(38, "cm"),
        barheight = grid::unit(0.65, "cm")
      )
    ) +
    scale_x_discrete(drop = TRUE, expand = expansion(add = 0)) +
    scale_y_continuous(
      breaks = meta$y_pos, labels = meta$sample_short,
      limits = c(0.5, nrow(meta) + 0.5), expand = expansion(mult = 0)
    ) +
    labs(x = NULL, y = NULL, fill = metric_axis_label(heatmap_col)) +
    heatmap_theme() +
    theme(
      strip.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13.8, face = "bold"),
      plot.margin = margin(5, 5, 5, 0)
    )

  if (!requireNamespace("patchwork", quietly = TRUE)) return(invisible(NULL))
  combined <- (dend | ann | tai_ann | hp) +
    patchwork::plot_layout(widths = c(2.1, 0.35, 0.48, 13), guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(
        "Hierarchical clustering of codon-level tRNA weights, annotated by ",
        clean_label(annotation_name)
      ),
      subtitle = NULL
    ) & theme(legend.position = "bottom", legend.box = "vertical")

  stem <- paste0("trna_weights_dendrogram_heatmap_by_", annotation_name)
  save_plot_pair(combined, stem, outdir, 30.0, max(heatmap_large_height, 6.5 + 0.26 * nrow(meta)), formats)
  readr::write_tsv(
    meta %>% select(sample, sample_short, annotation, cluster_order = y_pos),
    file.path(outdir, paste0(stem, "_sample_order.tsv"))
  )
}

readr::write_tsv(
  tibble(
    clustering_input = heatmap_col,
    distance = "1 - pairwise Spearman correlation between genome codon-weight profiles",
    linkage = "average",
    missing_value_handling = "codon-wise median imputation for clustering only",
    interpretation = "exploratory similarity dendrogram; not a phylogenetic reconstruction"
  ),
  file.path(outdir, "trna_weights_dendrogram_method.tsv")
)

make_dendrogram_heatmap(phylum)
make_dendrogram_heatmap(lifestyle)
readr::write_tsv(
  cluster_association_results,
  file.path(outdir, "trna_weights_annotation_association.tsv"),
  na = "NA"
)

message("Done: ", outdir)
