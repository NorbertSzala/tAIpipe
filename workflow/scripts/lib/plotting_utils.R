`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

# Shared static plotting utilities for ggplot2-based scripts.

read_plot_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read plotting config.", call. = FALSE)
  }
  yaml::read_yaml(path)
}

get_nested <- function(x, path, default = NULL) {
  current <- x
  for (key in path) {
    if (is.null(current[[key]])) {
      return(default)
    }
    current <- current[[key]]
  }
  current
}

project_theme <- function(cfg) {
  base_family <- get_nested(cfg, c("theme", "base_family"), default = "sans")
  base_size <- get_nested(cfg, c("theme", "base_size"), default = 10)

  ggplot2::theme_minimal(base_family = base_family, base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.position = get_nested(cfg, c("theme", "legend_position"), default = "right"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "plain")
    )
}

get_palette <- function(cfg, palette_name = NULL, n = NULL) {
  palettes <- get_nested(cfg, c("palettes"), default = list())

  if (is.null(palette_name)) {
    palette_name <- get_nested(cfg, c("defaults", "discrete_palette"), default = "okabe_ito_8")
  }

  colours <- palettes[[palette_name]]

  if (is.null(colours)) {
    stop("Palette not found in plotting config: ", palette_name, call. = FALSE)
  }

  colours <- unlist(colours, use.names = FALSE)

  if (!is.null(n) && n > length(colours)) {
    warning(
      "Requested ", n, " colours from palette '", palette_name,
      "', but only ", length(colours), " are defined. Colours will be reused.",
      call. = FALSE
    )
    colours <- rep(colours, length.out = n)
  }

  if (!is.null(n)) {
    colours <- colours[seq_len(n)]
  }

  colours
}

get_category_mapping <- function(cfg, mapping) {
  maps <- get_nested(cfg, c("category_mappings"), default = list())
  out <- maps[[mapping]]
  if (is.null(out)) {
    stop("Category mapping not found in plotting config: ", mapping, call. = FALSE)
  }
  unlist(out, use.names = TRUE)
}

scale_fill_project <- function(cfg, mapping = NULL, palette = NULL, levels = NULL, drop = FALSE, ...) {
  if (!is.null(mapping)) {
    values <- get_category_mapping(cfg, mapping)
  } else {
    n <- if (!is.null(levels)) length(levels) else NULL
    values <- get_palette(cfg, palette, n = n)
    if (!is.null(levels)) {
      names(values) <- levels
    }
  }

  ggplot2::scale_fill_manual(values = values, drop = drop, ...)
}

scale_colour_project <- function(cfg, mapping = NULL, palette = NULL, levels = NULL, drop = FALSE, ...) {
  if (!is.null(mapping)) {
    values <- get_category_mapping(cfg, mapping)
  } else {
    n <- if (!is.null(levels)) length(levels) else NULL
    values <- get_palette(cfg, palette, n = n)
    if (!is.null(levels)) {
      names(values) <- levels
    }
  }

  ggplot2::scale_colour_manual(values = values, drop = drop, ...)
}

scale_continuous_project <- function(cfg, option = "viridis", direction = 1, ...) {
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for continuous colour scales.", call. = FALSE)
  }
  viridis::scale_colour_viridis(option = option, direction = direction, ...)
}

scale_fill_continuous_project <- function(cfg, option = "viridis", direction = 1, ...) {
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for continuous fill scales.", call. = FALSE)
  }
  viridis::scale_fill_viridis(option = option, direction = direction, ...)
}

get_size_profile <- function(cfg, size = "wide") {
  profiles <- get_nested(cfg, c("dimensions"), default = list())
  profile <- profiles[[size]]

  # Backward-compatible defaults. Older config/plotting.yaml files did not
  # define dimensions even though several plotting scripts request standard
  # aliases such as "wide", "square" and "tall". Failing here stops otherwise
  # valid plotting jobs, so use conservative publication-size defaults and keep
  # config-provided values authoritative when present.
  if (is.null(profile)) {
    defaults <- list(
      single = list(width = 3.35, height = 2.70, units = "in"),
      double = list(width = 7.10, height = 4.80, units = "in"),
      wide = list(width = 10.00, height = 5.50, units = "in"),
      square = list(width = 7.10, height = 7.10, units = "in"),
      tall = list(width = 7.10, height = 9.00, units = "in"),
      tall_facets = list(width = 8.00, height = 12.00, units = "in")
    )
    profile <- defaults[[size]]
  }

  if (is.null(profile)) {
    stop("Size profile not found in plotting config and no fallback exists: ", size, call. = FALSE)
  }

  profile
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

# Apply the requested two-point increase after layout so that every visible
# text grob grows equally, regardless of the base theme used by a script.
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

save_plot <- function(plot, filename, cfg, size = "wide", dpi = NULL, device = NULL) {
  profile <- get_size_profile(cfg, size)

  width <- as.numeric(profile[["width"]])
  height <- as.numeric(profile[["height"]])
  units <- profile[["units"]] %||% get_nested(cfg, c("output", "units"), default = "in")
  if (is.null(dpi)) {
    dpi <- get_nested(cfg, c("output", "raster", "dpi"), default = 300)
  }

  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  ext <- tolower(tools::file_ext(filename))

  if (is.null(device)) {
    device <- switch(
      ext,
      png = "png",
      pdf = grDevices::cairo_pdf,
      svg = "svg",
      NULL
    )
  }

  ggplot2::ggsave(
    filename = filename,
    plot = prepare_plot_for_export(plot),
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    device = device,
    limitsize = FALSE
  )

  invisible(filename)
}

save_plot_outputs <- function(plot, outputs, cfg, size = "wide", dpi = NULL) {
  for (output in outputs) {
    save_plot(plot = plot, filename = output, cfg = cfg, size = size, dpi = dpi)
  }
  invisible(outputs)
}

combine_plots <- function(plots, cfg = NULL, layout = c("paired", "stacked", "grid_2x2", "grid_3col")) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required to combine plots.", call. = FALSE)
  }

  layout <- match.arg(layout)

  ncol <- switch(
    layout,
    paired = 2,
    stacked = 1,
    grid_2x2 = 2,
    grid_3col = 3
  )

  patchwork::wrap_plots(plots, ncol = ncol)
}

recommended_facet_size <- function(cfg, n_panels, ncol = NULL, panel_width = 3.2, panel_height = 2.6) {
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_panels))
  }

  nrow <- ceiling(n_panels / ncol)

  list(
    width = max(6, ncol * panel_width),
    height = max(4, nrow * panel_height),
    units = "in"
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


format_p_value <- function(x) {
    dplyr::case_when(
        is.na(x) ~ "p = NA",
        x < 0.001 ~ paste0("p = ", format(x, scientific = TRUE, digits = 2)),
        TRUE ~ paste0("p = ", formatC(x, format = "f", digits = 3))
    )
}
