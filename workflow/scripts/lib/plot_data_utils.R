# Data preparation helpers for plot-ready summary tables.

read_tsv_checked <- function(path, required_columns = NULL, table_name = path) {
  data <- readr::read_tsv(path, show_col_types = FALSE)
  require_non_empty(data, table_name)
  if (!is.null(required_columns)) {
    require_columns(data, required_columns, table_name)
  }
  data
}

write_tsv_safe <- function(data, path) {
  safe_make_dir(dirname(path))
  readr::write_tsv(data, path)
  invisible(path)
}

order_by_median <- function(data, group_col, value_col, decreasing = FALSE) {
  levels <- data |>
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[value_col]])) |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::summarise(.median = stats::median(.data[[value_col]], na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(if (decreasing) dplyr::desc(.median) else .median) |>
    dplyr::pull(.data[[group_col]])

  factor(data[[group_col]], levels = levels)
}

order_by_frequency <- function(data, group_col, decreasing = TRUE) {
  levels <- data |>
    dplyr::count(.data[[group_col]], name = ".n") |>
    dplyr::arrange(if (decreasing) dplyr::desc(.n) else .n) |>
    dplyr::pull(.data[[group_col]])

  factor(data[[group_col]], levels = levels)
}

select_top_n_by <- function(data, column, n = 20, decreasing = TRUE, with_ties = FALSE) {
  if (!column %in% names(data)) {
    stop("Column not found: ", column, call. = FALSE)
  }

  data |>
    dplyr::arrange(if (decreasing) dplyr::desc(.data[[column]]) else .data[[column]]) |>
    dplyr::slice_head(n = n)
}

select_top_terms <- function(data, p_col = "padj", n = 20) {
  if (!p_col %in% names(data)) {
    stop("p-value column not found: ", p_col, call. = FALSE)
  }

  data |>
    dplyr::filter(!is.na(.data[[p_col]])) |>
    dplyr::arrange(.data[[p_col]]) |>
    dplyr::slice_head(n = n)
}

summarize_binary_by_group <- function(data, group_cols, feature_map) {
  purrr::imap_dfr(
    feature_map,
    function(column_name, feature_label) {
      require_columns(data, c(group_cols, column_name), "binary feature input")

      data |>
        dplyr::mutate(.present = coerce_binary_feature(.data[[column_name]], column_name)) |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
        dplyr::summarise(
          feature = feature_label,
          n_total = sum(!is.na(.present)),
          n_present = sum(.present, na.rm = TRUE),
          proportion = dplyr::if_else(n_total > 0, n_present / n_total, NA_real_),
          .groups = "drop"
        )
    }
  )
}

summarize_numeric_by_group <- function(data, group_cols, feature_map) {
  purrr::imap_dfr(
    feature_map,
    function(column_name, feature_label) {
      require_columns(data, c(group_cols, column_name), "numeric feature input")
      check_numeric_columns(data, column_name, "numeric feature input")

      data |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
        dplyr::summarise(
          feature = feature_label,
          n = sum(!is.na(.data[[column_name]])),
          mean = mean(.data[[column_name]], na.rm = TRUE),
          median = stats::median(.data[[column_name]], na.rm = TRUE),
          sd = stats::sd(.data[[column_name]], na.rm = TRUE),
          q1 = as.numeric(stats::quantile(.data[[column_name]], 0.25, na.rm = TRUE)),
          q3 = as.numeric(stats::quantile(.data[[column_name]], 0.75, na.rm = TRUE)),
          .groups = "drop"
        )
    }
  )
}

collapse_rare_categories <- function(data, column, max_categories = 12, other_label = "Other") {
  require_columns(data, column, "category input")

  top_levels <- data |>
    dplyr::count(.data[[column]], sort = TRUE) |>
    dplyr::slice_head(n = max_categories) |>
    dplyr::pull(.data[[column]])

  data |>
    dplyr::mutate(
      "{column}" := dplyr::if_else(
        .data[[column]] %in% top_levels,
        as.character(.data[[column]]),
        other_label
      )
    )
}

warn_if_many_categories <- function(x, limit = 8, context = "plot") {
  n <- length(unique(stats::na.omit(x)))
  if (n > limit) {
    warning(
      context, " has ", n, " categories. Consider faceting, grouping rare categories, ",
      "or using shape/linetype in addition to colour.",
      call. = FALSE
    )
  }
  invisible(n)
}
