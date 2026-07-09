# Utility functions for validating tables used by statistics and plotting scripts.
# These functions should have no project-specific dependencies except base R.

require_columns <- function(data, columns, table_name = "table") {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(
      "Missing columns in ", table_name, ": ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

require_non_empty <- function(data, table_name = "table") {
  if (!is.data.frame(data)) {
    stop(table_name, " is not a data.frame.", call. = FALSE)
  }
  if (nrow(data) == 0) {
    stop(table_name, " is empty.", call. = FALSE)
  }
  invisible(TRUE)
}

check_numeric_columns <- function(data, columns, table_name = "table") {
  require_columns(data, columns, table_name)

  non_numeric <- columns[!vapply(data[columns], is.numeric, logical(1))]

  if (length(non_numeric) > 0) {
    stop(
      "Non-numeric columns in ", table_name, ": ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

check_no_missing_keys <- function(data, keys, table_name = "table") {
  require_columns(data, keys, table_name)

  bad <- keys[vapply(data[keys], function(x) any(is.na(x) | x == ""), logical(1))]

  if (length(bad) > 0) {
    stop(
      "Missing key values in ", table_name, ": ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

check_unique_keys <- function(data, keys, table_name = "table") {
  require_columns(data, keys, table_name)

  duplicated_rows <- duplicated(data[keys])

  if (any(duplicated_rows)) {
    example <- utils::head(data[duplicated_rows, keys, drop = FALSE], 5)
    stop(
      "Duplicated key rows in ", table_name, " for keys: ",
      paste(keys, collapse = ", "),
      "\nExamples:\n",
      paste(utils::capture.output(print(example)), collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

coerce_binary_feature <- function(x, column_name = "binary feature") {
  if (is.logical(x)) {
    return(x)
  }

  if (is.numeric(x) || is.integer(x)) {
    unique_values <- unique(stats::na.omit(x))
    if (!all(unique_values %in% c(0, 1))) {
      warning(
        "Column '", column_name,
        "' is numeric but not strictly 0/1. Values > 0 will be treated as TRUE.",
        call. = FALSE
      )
    }
    return(ifelse(is.na(x), NA, x > 0))
  }

  if (is.character(x) || is.factor(x)) {
    y <- tolower(trimws(as.character(x)))
    true_values <- c("true", "t", "1", "yes", "y", "present", "presence")
    false_values <- c("false", "f", "0", "no", "n", "absent", "absence")

    out <- rep(NA, length(y))
    out[y %in% true_values] <- TRUE
    out[y %in% false_values] <- FALSE

    unknown <- unique(y[!is.na(y) & !(y %in% c(true_values, false_values))])
    if (length(unknown) > 0) {
      stop(
        "Column '", column_name,
        "' cannot be coerced to binary. Unknown values: ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }

    return(out)
  }

  stop("Column '", column_name, "' is not binary-compatible.", call. = FALSE)
}

safe_make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}
