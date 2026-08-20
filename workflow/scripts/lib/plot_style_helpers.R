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
    mean_delta_ENC = "Mean delta ENC",
    median_delta_ENC = "Median delta ENC",
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

# Convert mixed plain text and inline LaTeX mathematics through latex2exp.
# Plain prose should remain outside $...$; only mathematical fragments belong
# in math mode. This avoids latex2exp parser failures caused by long
# \textrm{} blocks, reserved R words such as "in", and unsupported "\\ "
# spacing commands inside \mathrm{}.
plain_label_from_tex <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\$", "")
  x <- stringr::str_replace_all(x, "\\\\(?:textrm|text|mathrm|operatorname)\\{([^{}]*)\\}", "\\1")
  x <- stringr::str_replace_all(x, "\\\\overline\\{([^{}]*)\\}", "mean \\1")
  x <- stringr::str_replace_all(x, "\\\\Delta", "Delta ")
  x <- stringr::str_replace_all(x, "\\\\rho", "rho")
  x <- stringr::str_replace_all(x, "\\\\tau", "tau")
  x <- stringr::str_replace_all(x, "\\\\log_\\{([0-9]+)\\}", "log\\1")
  x <- stringr::str_replace_all(x, "\\\\(?:;|,|quad|qquad)", " ")
  x <- stringr::str_replace_all(x, "[{}]", "")
  stringr::str_squish(x)
}

tex_label <- function(x, output = "expression") {
  x <- as.character(x)
  if (!output %in% c("expression", "character")) {
    stop("Unsupported tex_label output: ", output)
  }

  if (requireNamespace("latex2exp", quietly = TRUE)) {
    rendered <- tryCatch(
      latex2exp::TeX(x, output = output),
      error = function(e) {
        warning(
          "latex2exp could not render label; using a plain-text fallback. ",
          "Label: ", paste(x, collapse = " | "), ". Error: ",
          conditionMessage(e),
          call. = FALSE,
          immediate. = TRUE
        )
        NULL
      }
    )
    if (!is.null(rendered)) return(rendered)
  }

  fallback <- plain_label_from_tex(x)
  if (identical(output, "character")) {
    return(vapply(
      fallback,
      function(value) encodeString(value, quote = "\""),
      character(1)
    ))
  }
  as.expression(fallback)
}

metric_tex_string <- function(metric, mean_symbol = FALSE) {
  metric <- as.character(metric)
  if (length(metric) != 1L) {
    return(vapply(metric, metric_tex_string, character(1), mean_symbol = mean_symbol))
  }
  metric <- metric[[1]]
  labels <- c(
    tAI = "$\\mathrm{tAI}$",
    tAI_z = "$\\mathrm{tAI}_{z}$",
    CAI = "$\\mathrm{CAI}$",
    ENC = "$\\mathrm{ENC}$",
    delta_ENC = "$\\Delta\\mathrm{ENC}$",
    GC = "$\\mathrm{GC}$",
    mean_GC = "$\\overline{\\mathrm{GC}}$",
    median_GC = "$\\mathrm{median}(\\mathrm{GC})$",
    genome_gc = "$\\mathrm{GC}_{\\mathrm{genome}}$",
    GC3s = "$\\mathrm{GC}_{3s}$",
    genome_RSCU = "$\\mathrm{RSCU}_{\\mathrm{genome}}$",
    RSCU = "$\\mathrm{RSCU}$",
    rscu = "$\\mathrm{RSCU}$",
    tRNA_weight = "$w_{\\mathrm{tRNA}}$",
    trna_weight = "$w_{\\mathrm{tRNA}}$",
    trna_absolute_weight = "$w_{\\mathrm{tRNA}}$",
    trna_copy_number = "$n_{\\mathrm{tRNA}}$",
    relative_trna = "Relative $\\mathrm{tRNA}$ abundance",
    codon_frequency = "Codon frequency",
    frequency = "Codon frequency",
    codon_count = "Codon count",
    mean_tAI = "$\\overline{\\mathrm{tAI}}$",
    median_tAI = "$\\mathrm{median}(\\mathrm{tAI})$",
    mean_CAI = "$\\overline{\\mathrm{CAI}}$",
    median_CAI = "$\\mathrm{median}(\\mathrm{CAI})$",
    mean_ENC = "$\\overline{\\mathrm{ENC}}$",
    median_ENC = "$\\mathrm{median}(\\mathrm{ENC})$",
    mean_delta_ENC = "$\\overline{\\Delta\\mathrm{ENC}}$",
    median_delta_ENC = "$\\mathrm{median}(\\Delta\\mathrm{ENC})$",
    mean_GC3s = "$\\overline{\\mathrm{GC}_{3s}}$",
    median_GC3s = "$\\mathrm{median}(\\mathrm{GC}_{3s})$"
  )
  if (metric %in% names(labels)) return(unname(labels[[metric]]))
  clean_label(metric)
}

# Device-independent labels used when plotmath/latex2exp is not worth the risk of
# literal TeX commands appearing in raster output. These labels preserve the
# biological meaning of means, medians and standardized quantities.
metric_plain_label <- function(metric) {
  metric <- as.character(metric)
  labels <- c(
    tAI = "tAI",
    tAI_z = "tAI z-score",
    CAI = "CAI",
    ENC = "ENC",
    delta_ENC = "Delta ENC",
    GC = "GC",
    mean_GC = "Mean GC",
    median_GC = "Median GC",
    genome_gc = "Genome GC",
    GC3s = "GC3s",
    genome_RSCU = "Genomic RSCU",
    RSCU = "RSCU",
    rscu = "RSCU",
    tRNA_weight = "tRNA weight",
    trna_weight = "tRNA weight",
    trna_absolute_weight = "Absolute tRNA weight",
    trna_copy_number = "tRNA copy number",
    relative_trna = "Relative tRNA abundance",
    codon_frequency = "Codon frequency",
    frequency = "Codon frequency",
    codon_count = "Codon count",
    mean_tAI = "Mean tAI",
    median_tAI = "Median tAI",
    mean_CAI = "Mean CAI",
    median_CAI = "Median CAI",
    mean_ENC = "Mean ENC",
    median_ENC = "Median ENC",
    mean_delta_ENC = "Mean Delta ENC",
    median_delta_ENC = "Median Delta ENC",
    mean_GC3s = "Mean GC3s",
    median_GC3s = "Median GC3s"
  )
  out <- unname(labels[metric])
  missing <- is.na(out)
  out[missing] <- clean_label(metric[missing])
  out
}

metric_axis_label <- function(metric) {
  metric_plain_label(metric)
}

metric_title_label <- function(metric) {
  metric_axis_label(metric)
}

metric_parsed_label <- function(metric) {
  metric_plain_label(metric)
}

metric_cv_label <- function(metric) {
  # Use a descriptive plain label for genomic RSCU and native plotmath for the
  # remaining compact CV labels. Neither route passes nested TeX through a
  # labeller, which can expose literal \mathrm text in raster output.
  native <- switch(as.character(metric),
    genome_RSCU = "Codon variation based on genomic RSCU",
    RSCU = expression(CV(RSCU)),
    rscu = expression(CV(RSCU)),
    tRNA_weight = expression(CV(w[tRNA])),
    trna_weight = expression(CV(w[tRNA])),
    trna_absolute_weight = expression(CV(W[tRNA])),
    trna_copy_number = expression(CV(n[tRNA])),
    NULL
  )
  if (!is.null(native)) return(native)

  metric_string <- metric_tex_string(metric)
  if (grepl("^\\$.*\\$$", metric_string)) {
    inner <- sub("^\\$", "", sub("\\$$", "", metric_string))
    return(tex_label(paste0("$\\mathrm{CV}(", inner, ")$")))
  }
  tex_label(paste0("CV of ", metric_string))
}

math_labels <- function(key) {
  values <- c(
    spearman_rho = "Spearman rho",
    pearson_r = "Pearson r",
    kendall_tau = "Kendall tau",
    log2_or_cmh = "log2(CMH odds ratio)",
    log2_or_cmh_ci = "log2(CMH odds ratio) [95% CI]",
    minus_log10_fdr = "-log10(FDR)",
    r_squared = "R-squared",
    z_score = "Within-genome z-score"
  )
  if (!key %in% names(values)) stop("Unknown mathematical label: ", key)
  unname(values[[key]])
}

correlation_tex_string <- function(variable) {
  variable <- as.character(variable)
  vapply(variable, function(v) {
    if (v %in% c(
      "tAI", "tAI_z", "CAI", "ENC", "delta_ENC", "GC", "mean_GC",
      "median_GC", "GC3s", "genome_gc", "mean_tAI", "median_tAI",
      "mean_CAI", "median_CAI", "mean_ENC", "median_ENC",
      "mean_delta_ENC", "median_delta_ENC", "mean_GC3s", "median_GC3s",
      "genome_RSCU", "RSCU", "rscu"
    )) return(metric_tex_string(v))
    switch(v,
      protein_length_aa = "Protein length $[\\mathrm{aa}]$",
      log_protein_length_aa = "$\\log(1+L_{\\mathrm{protein}})$",
      lcr_total_length = "Total LCR length $[\\mathrm{aa}]$",
      tm_total_length = "Total TM length $[\\mathrm{aa}]$",
      lcr_count = "LCR count",
      tm_count = "TM count",
      fraction_lcr = "LCR fraction",
      fraction_tm = "TM fraction",
      fraction_signal_peptide = "Signal-peptide fraction",
      fraction_pfam = "PFAM fraction",
      clean_label(v)
    )
  }, character(1))
}

correlation_axis_labels <- function(variable) {
  # Build native R plotmath expressions directly. Unlike escaped TeX strings,
  # these are parsed once and remain stable in PNG, PDF and SVG devices.
  plotmath_string <- function(v) {
    switch(as.character(v),
      tAI = 'plain("tAI")',
      tAI_z = 'plain("tAI")[z]',
      CAI = 'plain("CAI")',
      ENC = 'plain("ENC")',
      delta_ENC = 'Delta*plain("ENC")',
      GC = 'plain("GC")',
      mean_GC = 'bar(plain("GC"))',
      median_GC = 'median(plain("GC"))',
      genome_gc = 'plain("GC")[genome]',
      GC3s = 'plain("GC")[3*s]',
      mean_tAI = 'bar(plain("tAI"))',
      median_tAI = 'median(plain("tAI"))',
      mean_CAI = 'bar(plain("CAI"))',
      median_CAI = 'median(plain("CAI"))',
      mean_ENC = 'bar(plain("ENC"))',
      median_ENC = 'median(plain("ENC"))',
      mean_delta_ENC = 'bar(Delta*plain("ENC"))',
      median_delta_ENC = 'median(Delta*plain("ENC"))',
      mean_GC3s = 'bar(plain("GC")[3*s])',
      median_GC3s = 'median(plain("GC")[3*s])',
      protein_length_aa = 'plain("Protein length [aa]")',
      log_protein_length_aa = 'log(1+L[protein])',
      lcr_total_length = 'plain("Total LCR length [aa]")',
      tm_total_length = 'plain("Total TM length [aa]")',
      lcr_count = 'plain("LCR count")',
      tm_count = 'plain("TM count")',
      fraction_lcr = 'plain("Genes with LCR [%]")',
      fraction_tm = 'plain("Genes with TM [%]")',
      fraction_signal_peptide = 'plain("Genes with signal peptide [%]")',
      fraction_pfam = 'plain("Genes with PFAM [%]")',
      paste0("plain(", encodeString(clean_label(v), quote = '"'), ")")
    )
  }
  as.expression(lapply(vapply(variable, plotmath_string, character(1)), function(x) {
    parse(text = x, keep.source = FALSE)[[1]]
  }))
}

correlation_plain_labels <- function(variable) {
  variable <- as.character(variable)
  vapply(variable, function(v) {
    if (v %in% c(
      "tAI", "tAI_z", "CAI", "ENC", "delta_ENC", "GC", "mean_GC",
      "median_GC", "GC3s", "genome_gc", "mean_tAI", "median_tAI",
      "mean_CAI", "median_CAI", "mean_ENC", "median_ENC",
      "mean_delta_ENC", "median_delta_ENC", "mean_GC3s", "median_GC3s",
      "genome_RSCU", "RSCU", "rscu"
    )) return(metric_plain_label(v))
    switch(v,
      protein_length_aa = "Protein length [aa]",
      log_protein_length_aa = "log(1 + protein length)",
      lcr_total_length = "Total LCR length [aa]",
      tm_total_length = "Total TM length [aa]",
      lcr_count = "LCR count",
      tm_count = "TM count",
      fraction_lcr = "Genes with LCR [%]",
      fraction_tm = "Genes with TM [%]",
      fraction_signal_peptide = "Genes with signal peptide [%]",
      fraction_pfam = "Genes with PFAM [%]",
      clean_label(v)
    )
  }, character(1))
}

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
  # Vectorised formatter: accepts one p-value or an entire numeric column.
  # Exact p-values remain available in the corresponding TSV result tables.
  p_num <- suppressWarnings(as.numeric(as.character(p)))

  out <- rep("p = NA", length(p_num))

  finite <- is.finite(p_num)

  out[finite & p_num < 0.001] <- "p < 0.001"
  out[finite & p_num >= 0.001] <- paste0(
    "p = ", formatC(p_num[finite & p_num >= 0.001], format = "f", digits = 3)
  )

  names(out) <- names(p)
  out
}

sig_label <- function(p) {
  p_label(p)
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

# Paired two-level comparison after reshaping by an explicit observational unit
# such as genome. `level_order` is c(reference, comparison).
safe_paired_wilcox_p <- function(data, id, group, value, level_order) {
  id_q <- rlang::enquo(id)
  group_q <- rlang::enquo(group)
  value_q <- rlang::enquo(value)
  if (length(level_order) != 2L) return(NA_real_)

  wide <- data %>%
    dplyr::transmute(
      .id = as.character(!!id_q),
      .group = as.character(!!group_q),
      .value = suppressWarnings(as.numeric(!!value_q))
    ) %>%
    dplyr::filter(!is.na(.id), nzchar(.id), .group %in% level_order, is.finite(.value)) %>%
    dplyr::group_by(.id, .group) %>%
    dplyr::summarise(.value = stats::median(.value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = .group, values_from = .value)

  if (!all(level_order %in% names(wide))) return(NA_real_)
  ok <- is.finite(wide[[level_order[[1]]]]) & is.finite(wide[[level_order[[2]]]])
  if (sum(ok) < 4L) return(NA_real_)
  tryCatch(
    stats::wilcox.test(
      wide[[level_order[[2]]]][ok], wide[[level_order[[1]]]][ok],
      paired = TRUE, exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

lm_annotation <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 5L || length(unique(x)) < 3L) {
    return("y = NA\nR-squared = NA\np = NA")
  }

  fit <- tryCatch(stats::lm(y ~ x), error = function(e) NULL)
  if (is.null(fit)) {
    return("y = NA\nR-squared = NA\np = NA")
  }

  sm <- summary(fit)
  intercept <- unname(stats::coef(fit)[[1]])
  slope <- unname(stats::coef(fit)[[2]])
  p <- tryCatch(sm$coefficients[2, 4], error = function(e) NA_real_)
  sign_text <- if (slope < 0) " - " else " + "
  equation_text <- paste0(
    "y = ", formatC(intercept, format = "f", digits = 3),
    sign_text, formatC(abs(slope), format = "f", digits = 3), "x"
  )
  r2_text <- paste0(
    "R-squared = ", formatC(sm$r.squared, format = "f", digits = 3)
  )
  p_text <- if (!is.finite(p)) {
    "p = NA"
  } else if (p < 1e-6) {
    "p < 1e-6"
  } else {
    paste0("p = ", formatC(p, format = "g", digits = 3))
  }
  paste(equation_text, r2_text, p_text, sep = "\n")
}

annotate_top_right <- function(label, size = 3.0) {
  parse_label <- inherits(label, "plotmath_label")
  ggplot2::annotate(
    "text",
    x = Inf,
    y = Inf,
    label = as.character(label),
    hjust = 1.04,
    vjust = 1.08,
    size = size,
    lineheight = 0.95,
    parse = parse_label
  )
}

remove_generated_plot_annotations <- function(plot) {
  if (inherits(plot, "patchwork")) {
    plot <- plot + patchwork::plot_annotation(
      title = NULL,
      subtitle = NULL,
      caption = NULL,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_blank(),
        plot.caption = ggplot2::element_blank(),
        plot.tag = ggplot2::element_blank()
      )
    )
    return(
      plot & ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_blank(),
        plot.caption = ggplot2::element_blank(),
        plot.tag = ggplot2::element_blank()
      )
    )
  }
  if (inherits(plot, "ggplot")) {
    return(
      plot +
        ggplot2::labs(
          title = NULL,
          subtitle = NULL,
          caption = NULL,
          tag = NULL,
          alt = NULL,
          alt_insight = NULL
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_blank(),
          plot.subtitle = ggplot2::element_blank(),
          plot.caption = ggplot2::element_blank(),
          plot.tag = ggplot2::element_blank()
        )
    )
  }
  plot
}

# Increase text only after a plot has been built. This preserves the relative
# sizes specified by individual plotting scripts while adding two points to all
# visible text, including axes, legends, facet strips and statistical labels.
increase_grob_font_size <- function(grob, points = 2) {
  if (inherits(grob, "text")) {
    font_size <- grob$gp$fontsize
    if (!is.null(font_size) && all(is.finite(font_size))) {
      grob$gp$fontsize <- font_size + points
    }
  }

  if (!is.null(grob$children)) {
    for (i in seq_along(grob$children)) {
      grob$children[[i]] <- increase_grob_font_size(grob$children[[i]], points)
    }
  }
  if (!is.null(grob$grobs)) {
    for (i in seq_along(grob$grobs)) {
      grob$grobs[[i]] <- increase_grob_font_size(grob$grobs[[i]], points)
    }
  }
  grob
}

prepare_plot_for_export <- function(plot, font_increase = 2) {
  plot <- remove_generated_plot_annotations(plot)
  grob <- if (inherits(plot, "patchwork")) {
    patchwork::patchworkGrob(plot)
  } else if (inherits(plot, "ggplot")) {
    ggplot2::ggplotGrob(plot)
  } else {
    plot
  }
  increase_grob_font_size(grob, font_increase)
}

save_plot_pair <- function(plot, stem, output_dir, width, height, formats = c("png", "pdf"), dpi = 300) {
  plot <- prepare_plot_for_export(plot)
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
  # Produce human-readable organism labels for plot axes and sample-order TSVs.
  # Rules:
  #   Neurospora_crassa_OR74A -> N. crassa OR74A when a strain is needed
  #   Basidiobolus_helicus_helicus -> B. helicus (duplicated species epithet removed)
  #   Haplosporangium_sp_Z_27 -> Haplosporangium sp Z 27 (do not abbreviate genus before "sp")
  samples <- as.character(samples)
  parts <- strsplit(samples, "_", fixed = TRUE)
  genus <- vapply(parts, function(x) if (length(x) >= 1) x[[1]] else NA_character_, character(1))
  species <- vapply(parts, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))

  base <- ifelse(
    !is.na(genus) & !is.na(species) & species == "sp",
    paste(genus, "sp"),
    ifelse(!is.na(genus) & !is.na(species), paste0(substr(genus, 1, 1), ". ", species), samples)
  )
  dup_base <- base %in% base[duplicated(base)]

  strain <- vapply(seq_along(parts), function(i) {
    x <- parts[[i]]
    if (length(x) <= 2) return(NA_character_)
    tail <- x[3:length(x)]
    tail <- tail[tail != species[[i]]]
    if (length(tail) == 0L) NA_character_ else paste(tail, collapse = " ")
  }, character(1))

  # Always keep suffixes for "Genus sp ..." labels, because otherwise the clade
  # is clear but the taxon identity becomes ambiguous.
  keep_suffix <- (!is.na(species) & species == "sp") | dup_base
  out <- ifelse(keep_suffix & !is.na(strain) & nzchar(strain), paste(base, strain), base)
  out <- stringr::str_replace_all(out, "_", " ")
  out <- stringr::str_squish(out)
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

order_project_groups <- function(groups, group_variable = NULL) {
  groups_chr <- unique(as.character(groups))
  priority <- character()
  if (identical(as.character(group_variable), "phylum") || any(c("Ascomycota", "Basidiomycota") %in% groups_chr)) {
    priority <- c("Ascomycota", "Basidiomycota")
  }
  c(intersect(priority, groups_chr), sort(setdiff(groups_chr, priority)))
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


binary_grey_values <- function(labels = c("Absent", "Present")) {
  vals <- c("Absent" = "#B8B8B8", "Present" = "#1F4E79",
            "absence" = "#B8B8B8", "presence" = "#1F4E79",
            "FALSE" = "#B8B8B8", "TRUE" = "#1F4E79",
            "no PFAM-LCR overlap" = "#B8B8B8", "PFAM-LCR overlap" = "#1F4E79")
  out <- vals[labels]
  missing <- is.na(out)
  out[missing] <- rep(c("#B8B8B8", "#1F4E79"), length.out = sum(missing))
  names(out) <- labels
  out
}


# Split-violin geom without an additional package. The first fill level is drawn
# on the left and the second on the right. Use a constant x within each feature.
GeomSplitViolin <- ggplot2::ggproto(
  "GeomSplitViolin", ggplot2::GeomViolin,
  default_aes = ggplot2::aes(
    weight = 1, colour = "grey20", fill = "white", linewidth = 0.5,
    alpha = NA, linetype = "solid", split_side = NA_real_
  ),
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data$group[1]
    mapped_side <- suppressWarnings(as.numeric(data$split_side[1]))
    draw_left <- if (is.finite(mapped_side)) mapped_side < 0 else grp %% 2 == 1
    newdata <- transform(data, x = if (draw_left) xminv else xmaxv)
    newdata <- newdata[order(if (draw_left) newdata$y else -newdata$y), , drop = FALSE]
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0L && !scales::zero_range(range(data$y))) {
      stop("draw_quantiles is not supported by geom_split_violin_project()")
    }
    ggplot2::GeomPolygon$draw_panel(newdata, ...)
  }
)

split_side_from_level <- function(x, left_level) {
  ifelse(as.character(x) == as.character(left_level), -1, 1)
}

geom_split_violin_project <- function(mapping = NULL, data = NULL, stat = "ydensity",
                                      position = "identity", ..., trim = TRUE,
                                      scale = "area", na.rm = FALSE,
                                      show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, na.rm = na.rm, ...)
  )
}

# Combine a scatterplot with group-wise marginal densities using patchwork.
# This avoids a hard dependency on ggExtra and keeps group colours consistent.
add_marginal_densities <- function(main_plot, data, x, y, group, colours,
                                   x_title = NULL, y_title = NULL,
                                   top_height = 1.2, right_width = 1.2) {
  if (!requireNamespace("patchwork", quietly = TRUE)) return(main_plot)
  xq <- rlang::enquo(x); yq <- rlang::enquo(y); gq <- rlang::enquo(group)
  top <- ggplot2::ggplot(data, ggplot2::aes(x = !!xq, fill = !!gq, colour = !!gq)) +
    ggplot2::geom_density(alpha = 0.22, linewidth = 0.55, na.rm = TRUE, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = colours, drop = FALSE) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::guides(fill = "none", colour = "none") +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL, colour = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(0, 0, 0, 0))
  right <- ggplot2::ggplot(data, ggplot2::aes(x = !!yq, fill = !!gq, colour = !!gq)) +
    ggplot2::geom_density(alpha = 0.22, linewidth = 0.55, na.rm = TRUE, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = colours, drop = FALSE) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::guides(fill = "none", colour = "none") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL, colour = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(0, 0, 0, 0))
  spacer <- patchwork::plot_spacer()
  (top + spacer + main_plot + right) +
    patchwork::plot_layout(widths = c(8, right_width), heights = c(top_height, 8), guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

# Marginal densities for a single ungrouped scatterplot.
add_single_marginal_densities <- function(main_plot, data, x, y,
                                          colour = "#1F4E79",
                                          top_height = 1.2,
                                          right_width = 1.2) {
  if (!requireNamespace("patchwork", quietly = TRUE)) return(main_plot)
  xq <- rlang::enquo(x)
  yq <- rlang::enquo(y)
  top <- ggplot2::ggplot(data, ggplot2::aes(x = !!xq)) +
    ggplot2::geom_density(fill = colour, colour = colour, alpha = 0.25,
                          linewidth = 0.55, na.rm = TRUE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  right <- ggplot2::ggplot(data, ggplot2::aes(x = !!yq)) +
    ggplot2::geom_density(fill = colour, colour = colour, alpha = 0.25,
                          linewidth = 0.55, na.rm = TRUE) +
    ggplot2::coord_flip() +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  (top + patchwork::plot_spacer() + main_plot + right) +
    patchwork::plot_layout(
      widths = c(8, right_width),
      heights = c(top_height, 8)
    )
}

# Add compact boxplots centred in the left/right half of split violins.
add_split_violin_boxplots <- function(width = 0.075, dodge = 0.16,
                                      outlier.shape = NA,
                                      fill = "white", colour = "grey25",
                                      alpha = 0.92, linewidth = 0.35) {
  ggplot2::geom_boxplot(
    width = width,
    position = ggplot2::position_dodge(width = dodge),
    outlier.shape = outlier.shape,
    fill = fill,
    colour = colour,
    alpha = alpha,
    linewidth = linewidth
  )
}

# Correlation-variable order equivalent in spirit to ggcorrplot(hc.order=TRUE):
# distance is 1 - |rho| and complete-linkage clustering determines the axis order.
# The order is deterministic for a fixed matrix.
cluster_correlation_variables <- function(cor_mat) {
  cor_mat <- as.matrix(cor_mat)
  if (nrow(cor_mat) < 3L) return(rownames(cor_mat))
  cor_mat[!is.finite(cor_mat)] <- 0
  diag(cor_mat) <- 1
  dmat <- 1 - abs(cor_mat)
  dmat[!is.finite(dmat)] <- 1
  dmat[dmat < 0] <- 0
  dmat <- (dmat + t(dmat)) / 2
  diag(dmat) <- 0
  d <- stats::as.dist(dmat)
  hc <- stats::hclust(d, method = "complete")
  hc$labels[hc$order]
}
