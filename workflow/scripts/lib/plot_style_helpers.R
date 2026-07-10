# Shared plotting helpers for tAIpipe.
# Keep only small, reusable formatting/statistics utilities here.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

wrap_text <- function(x, width = 90) {
  stringr::str_wrap(as.character(x), width = width)
}

clean_label <- function(x) {
  x <- as.character(x)
  recode <- c(
    lcr_present = "LCR",
    pfam_present = "PFAM domains",
    tm_present = "Transmembrane region",
    signal_peptide_present = "Signal peptide",
    pfam_lcr_present = "PFAM-LCR overlap",
    pfam_lcr_overlap_present = "PFAM-LCR overlap",
    algal_parasite = "Algal parasite",
    animal_endosymbiont = "Animal endosymbiont",
    animal_parasite = "Animal parasite",
    arbuscular_mycorrhizal = "Arbuscular mycorrhizal",
    ectomycorrhizal = "Ectomycorrhizal",
    lichenized = "Lichenized",
    plant_pathogen = "Plant pathogen",
    soil_saprotroph = "Soil saprotroph",
    unspecified_saprotroph = "Unspecified saprotroph",
    wood_saprotroph = "Wood saprotroph",
    nectar_tap_saprotroph = "Nectar tap saprotroph",
    pollen_saprotroph = "Pollen saprotroph",
    litter_saprotroph = "Litter saprotroph",
    dung_saprotroph = "Dung saprotroph",
    mycoparasite = "Mycoparasite",
    protistan_parasite = "Protistan parasite",
    mean_tAI = "Mean tAI",
    median_tAI = "Median tAI",
    mean_GC3s = "Mean GC3s",
    median_GC3s = "Median GC3s",
    genome_gc = "Genome GC",
    mean_CAI = "Mean CAI",
    median_CAI = "Median CAI",
    mean_ENC = "Mean ENC",
    median_ENC = "Median ENC",
    phylum = "Phylum",
    lifestyle = "Lifestyle"
  )
  out <- recode[x]
  missing <- is.na(out)

  out[missing] <- x[missing] %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish() %>%
    stringr::str_to_sentence()
  unname(out)
}

# Backward-compatible alias used by older plotting scripts.
clean_plot_label <- clean_label

coerce_binary_plot <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x > 0)
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("1", "true", "t", "yes", "y", "present") ~ TRUE,
    x %in% c("0", "false", "f", "no", "n", "absent") ~ FALSE,
    TRUE ~ NA
  )
}

p_label <- function(p) {
  if (!is.finite(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p, format = "f", digits = 3))
}

sig_label <- function(p) {
  if (!is.finite(p)) return("test: NA")
  paste0(ifelse(p < 0.05, "significant", "not significant"), "\n", p_label(p))
}

safe_wilcox_p <- function(y, group) {
  ok <- is.finite(y) & !is.na(group)
  y <- y[ok]
  group <- droplevels(as.factor(group[ok]))
  if (length(unique(group)) != 2L || length(y) < 4L) return(NA_real_)
  tryCatch(stats::wilcox.test(y ~ group)$p.value, error = function(e) NA_real_)
}

safe_kruskal_p <- function(y, group) {
  ok <- is.finite(y) & !is.na(group)
  y <- y[ok]
  group <- droplevels(as.factor(group[ok]))
  if (length(unique(group)) < 2L || length(y) < 4L) return(NA_real_)
  tryCatch(stats::kruskal.test(y ~ group)$p.value, error = function(e) NA_real_)
}

lm_annotation <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 5L || length(unique(x)) < 3L) return("slope = NA\nR² = NA\np = NA")
  fit <- tryCatch(stats::lm(y ~ x), error = function(e) NULL)
  if (is.null(fit)) return("slope = NA\nR² = NA\np = NA")
  sm <- summary(fit)
  slope <- unname(stats::coef(fit)[[2]])
  p <- tryCatch(sm$coefficients[2, 4], error = function(e) NA_real_)
  paste0(
    "slope = ", formatC(slope, format = "g", digits = 3),
    "\nR² = ", formatC(sm$r.squared, format = "f", digits = 3),
    "\n", p_label(p)
  )
}

annotate_top_right <- function(label, size = 3.0) {
  ggplot2::annotate(
    "text", x = Inf, y = Inf, label = label,
    hjust = 1.04, vjust = 1.08, size = size, lineheight = 0.95
  )
}

save_plot_pair <- function(plot, stem, output_dir, width, height, formats = c("png", "pdf"), dpi = 300) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (fmt in formats) {
    path <- file.path(output_dir, paste0(stem, ".", fmt))
    ggplot2::ggsave(path, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  }
}

# Shortens sample labels to G. species, adding strain only when the same base
# species appears more than once. Repeated species epithets are not treated as strains:
#   Basidiobolus_helicus_helicus -> B. helicus
#   Ordospora_colligata_OC4      -> O. colligata_OC4 only if O. colligata has duplicates.
short_sample_labels <- function(samples) {
  samples <- as.character(samples)
  parts <- strsplit(samples, "_", fixed = TRUE)
  genus <- vapply(parts, function(x) if (length(x) >= 1) x[[1]] else NA_character_, character(1))
  species <- vapply(parts, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
  base <- ifelse(
    !is.na(genus) & !is.na(species),
    paste0(substr(genus, 1, 1), ". ", species),
    samples
  )
  dup_base <- base %in% base[duplicated(base)]
  strain <- vapply(seq_along(parts), function(i) {
    x <- parts[[i]]
    if (length(x) <= 2) return(NA_character_)
    tail <- x[3:length(x)]
    tail <- tail[tail != species[[i]]]
    if (length(tail) == 0L) NA_character_ else paste(tail, collapse = "_")
  }, character(1))
  out <- ifelse(dup_base & !is.na(strain) & nzchar(strain), paste0(base, "_", strain), base)
  names(out) <- samples
  out
}

add_group_n <- function(data, group_col) {
  tab <- data |>
    dplyr::count(.data[[group_col]], name = "n") |>
    dplyr::mutate(label = paste0(clean_label(.data[[group_col]]), " (", n, ")"))
  labels <- stats::setNames(tab$label, tab[[group_col]])
  list(data = data, labels = labels)
}

collapse_small_groups <- function(data, group_col, min_n = 5, other_label = "Other") {
  tab <- data |>
    dplyr::count(.data[[group_col]], name = "n")
  keep <- tab |>
    dplyr::filter(n >= min_n) |>
    dplyr::pull(.data[[group_col]])
  other_members <- tab |>
    dplyr::filter(!.data[[group_col]] %in% keep) |>
    dplyr::pull(.data[[group_col]])
  out <- data |>
    dplyr::mutate("{group_col}" := ifelse(.data[[group_col]] %in% keep, .data[[group_col]], other_label))
  attr(out, "other_members") <- other_members
  out
}

# Deterministic named colours. If you need absolute permanence across all scripts,
# use these named vectors, not unnamed palettes whose mapping changes with categories.
project_category_colours <- function(categories) {
  categories <- sort(unique(as.character(categories)))
  base <- c(
    "Ascomycota" = "#0072B2",
    "Basidiomycota" = "#D55E00",
    "Microsporidia" = "#009E73",
    "Mucoromycota" = "#CC79A7",
    "Mortierellomycota" = "#E69F00",
    "Glomeromycota" = "#56B4E9",
    "Chytridiomycota" = "#F0E442",
    "Blastocladiomycota" = "#999999",
    "Animal parasite" = "#D55E00",
    "Plant pathogen" = "#009E73",
    "Algal parasite" = "#56B4E9",
    "Animal endosymbiont" = "#CC79A7",
    "Unspecified saprotroph" = "#999999",
    "Soil saprotroph" = "#E69F00",
    "Wood saprotroph" = "#0072B2",
    "Litter saprotroph" = "#F0E442",
    "Other" = "#BDBDBD"
  )
  missing <- setdiff(categories, names(base))
  if (length(missing) > 0L) {
    palette <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#BBBBBB")
    base[missing] <- rep(palette, length.out = length(missing))
  }
  base[categories]
}
