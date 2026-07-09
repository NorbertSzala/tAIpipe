# Label and formatting helpers for static plots.

format_p_value <- function(x, digits = 3) {
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x < 0.001] <- "<0.001"
  out[!is.na(x) & x >= 0.001] <- formatC(x[!is.na(x) & x >= 0.001], format = "f", digits = digits)
  out
}

format_fdr <- function(x, digits = 3) {
  out <- format_p_value(x, digits = digits)
  paste0("FDR = ", out)
}

format_percent <- function(x, digits = 1, input = c("fraction", "percent")) {
  input <- match.arg(input)
  value <- if (input == "fraction") 100 * x else x
  paste0(formatC(value, format = "f", digits = digits), "%")
}

clean_feature_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\baa\\b", "AA", x, ignore.case = TRUE)
  x <- gsub("\\btai\\b", "tAI", x, ignore.case = TRUE)
  x <- gsub("\\bcai\\b", "CAI", x, ignore.case = TRUE)
  x <- gsub("\\bgc3s\\b", "GC3s", x, ignore.case = TRUE)
  x <- gsub("\\bpfam\\b", "PFAM", x, ignore.case = TRUE)
  x <- gsub("\\blcr\\b", "LCR", x, ignore.case = TRUE)
  x <- gsub("\\btm\\b", "TM", x, ignore.case = TRUE)
  trimws(gsub("\\s+", " ", x))
}

wrap_labels <- function(x, width = 35) {
  vapply(
    x,
    function(z) paste(strwrap(z, width = width), collapse = "\n"),
    character(1)
  )
}

make_tail_group_levels <- function(tail_fractions) {
  percentages <- format(100 * as.numeric(tail_fractions), trim = TRUE, scientific = FALSE)
  c(
    "All genes",
    as.vector(rbind(
      paste0("Bottom ", percentages, "%"),
      paste0("Top ", percentages, "%")
    ))
  )
}
