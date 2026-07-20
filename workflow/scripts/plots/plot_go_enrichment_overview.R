# Creates the canonical GO-enrichment dotplot and confidence-interval forest
# plot, using namespace facets and neutral term styling for readability.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

source(snakemake@input[["plotting_utils"]])
source("workflow/scripts/lib/table_validation_utils.R")
source("workflow/scripts/lib/plot_data_utils.R")
source("workflow/scripts/lib/label_utils.R")
source("workflow/scripts/lib/plot_style_helpers.R")

cfg <- read_plot_config(snakemake@input[["plotting_config"]])

pick_first_existing <- function(data, candidates, label, required = TRUE) {
  hit <- intersect(candidates, names(data))
  if (length(hit) == 0L) {
    if (!required) return(NA_character_)
    stop(label, " column not found. Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[[1]]
}

terms <- read_tsv_checked(snakemake@input[["top_terms"]], table_name = "go_enrichment_top_terms.tsv")
p_col <- pick_first_existing(terms, c("q_value", "padj", "fdr", "FDR", "adjusted_p_value", "p_adjust", "p_adj"), "Adjusted p-value/FDR")
go_col <- pick_first_existing(terms, c("go_id", "go_term", "go_terms"), "GO ID")
name_col <- pick_first_existing(terms, c("go_name", "term_name", "description", "name"), "GO description", required = FALSE)
namespace_col <- pick_first_existing(terms, c("go_namespace", "namespace", "ontology"), "GO namespace", required = FALSE)
tail_col <- pick_first_existing(terms, c("tail", "tai_tail", "direction"), "tAI tail", required = FALSE)

or_col <- pick_first_existing(terms, c("common_odds_ratio", "odds_ratio", "OR", "estimate"), "Odds ratio")
lo_col <- pick_first_existing(terms, c("conf_low", "ci_low", "lower_ci", "lower"), "Lower confidence limit", required = FALSE)
hi_col <- pick_first_existing(terms, c("conf_high", "ci_high", "upper_ci", "upper"), "Upper confidence limit", required = FALSE)

namespace_colours <- c(
  "Biological process" = "#252525",
  "Molecular function" = "#777777",
  "Cellular component" = "#C8C8C8"
)
tail_shapes <- c("High tAI" = 21, "Low tAI" = 24)
namespace_facet_labeller <- labeller(
  namespace_plot = function(x) stringr::str_wrap(x, width = 12)
)

terms <- terms %>%
  mutate(
    go_id_plot = toupper(str_squish(as.character(.data[[go_col]]))),
    go_name_plot = if (!is.na(name_col)) str_squish(as.character(.data[[name_col]])) else "Description unavailable",
    namespace_plot = if (!is.na(namespace_col)) clean_label(.data[[namespace_col]]) else NA_character_,
    tail_plot = if (!is.na(tail_col)) case_when(
      str_detect(tolower(as.character(.data[[tail_col]])), "high") ~ "High tAI",
      str_detect(tolower(as.character(.data[[tail_col]])), "low") ~ "Low tAI",
      TRUE ~ NA_character_
    ) else "GO enrichment",
    adjusted_value = suppressWarnings(as.numeric(.data[[p_col]])),
    odds_ratio = suppressWarnings(as.numeric(.data[[or_col]])),
    conf_low_or = if (!is.na(lo_col)) suppressWarnings(as.numeric(.data[[lo_col]])) else NA_real_,
    conf_high_or = if (!is.na(hi_col)) suppressWarnings(as.numeric(.data[[hi_col]])) else NA_real_,
    log2_or = if_else(is.finite(odds_ratio) & odds_ratio > 0, log2(odds_ratio), NA_real_),
    log2_ci_low = if_else(is.finite(conf_low_or) & conf_low_or > 0, log2(conf_low_or), NA_real_),
    log2_ci_high = if_else(is.finite(conf_high_or) & conf_high_or > 0, log2(conf_high_or), NA_real_),
    minus_log10_fdr = pmin(-log10(pmax(adjusted_value, .Machine$double.xmin)), 320),
    term_label_plain = paste0(go_id_plot, " — ", if_else(is.na(go_name_plot) | !nzchar(go_name_plot), "Description unavailable", go_name_plot))
  ) %>%
  filter(
    str_detect(go_id_plot, "^GO:[0-9]{7}$"),
    !is.na(namespace_plot), nzchar(namespace_plot),
    !tolower(namespace_plot) %in% c("na", "n/a", "none", "unknown", "not available"),
    is.finite(log2_or), is.finite(minus_log10_fdr)
  )

if (nrow(terms) == 0L) stop("No finite, labelled GO enrichment rows are available for plotting.", call. = FALSE)

term_levels <- terms %>% arrange(namespace_plot, tail_plot, log2_or) %>% pull(term_label_plain) %>% unique() %>% rev()
terms <- terms %>% mutate(term_label = factor(term_label_plain, levels = term_levels))

common_theme <- theme_minimal(base_size = 14) +
  theme(
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 13.5, lineheight = 1.0),
    strip.background = element_rect(fill = "grey94", colour = "grey60", linewidth = 0.35),
    strip.placement = "outside",
    axis.text.y = element_text(size = 12.5, lineheight = 1.02),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 13.5),
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 12.5, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.margin = margin(t = 8),
    plot.margin = margin(8, 10, 8, 8)
  )

dotplot <- ggplot(terms, aes(x = log2_or, y = term_label, size = minus_log10_fdr, shape = tail_plot)) +
  geom_vline(xintercept = 0, linewidth = 0.45, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.92, colour = "grey15", fill = "#1F4E79", stroke = 0.50) +
  facet_grid(
    namespace_plot ~ ., scales = "free_y", space = "free_y", switch = "y",
    labeller = namespace_facet_labeller
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 52)) +
  scale_shape_manual(values = tail_shapes, drop = FALSE) +
  scale_size_continuous(range = c(5.0, 14.0)) +
  labs(
    x = math_labels("log2_or_cmh"), y = NULL,
    shape = "tAI tail", size = math_labels("minus_log10_fdr"),
    title = "Top GO enrichment terms",
    subtitle = NULL
  ) + common_theme +
  guides(
    shape = guide_legend(order = 1, nrow = 1, override.aes = list(size = 5)),
    size = guide_legend(order = 2, nrow = 1)
  )

forest <- ggplot(terms, aes(x = log2_or, y = term_label, shape = tail_plot)) +
  geom_vline(xintercept = 0, linewidth = 0.45, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(
    data = terms %>% filter(is.finite(log2_ci_low), is.finite(log2_ci_high)),
    aes(xmin = log2_ci_low, xmax = log2_ci_high), height = 0.18, linewidth = 0.72,
    colour = "#1F4E79"
  ) +
  geom_point(size = 5.4, alpha = 0.94, colour = "#1F4E79") +
  facet_grid(
    namespace_plot ~ ., scales = "free_y", space = "free_y", switch = "y",
    labeller = namespace_facet_labeller
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 52)) +
  scale_shape_manual(values = c("High tAI" = 16, "Low tAI" = 17), drop = FALSE) +
  labs(
    x = math_labels("log2_or_cmh_ci"), y = NULL,
    shape = "tAI tail",
    title = "GO enrichment effect sizes",
    subtitle = NULL
  ) + common_theme +
  guides(
    shape = guide_legend(order = 1, nrow = 1, override.aes = list(size = 5.4))
  )

save_plot(dotplot, snakemake@output[["dotplot_png"]], cfg, size = "go_extra_wide")
save_plot(dotplot, snakemake@output[["dotplot_pdf"]], cfg, size = "go_extra_wide")
save_plot(forest, snakemake@output[["forest_png"]], cfg, size = "go_extra_wide")
save_plot(forest, snakemake@output[["forest_pdf"]], cfg, size = "go_extra_wide")
