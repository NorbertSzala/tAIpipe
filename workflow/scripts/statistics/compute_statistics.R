#!/usr/bin/env Rscript
# Legacy combined statistical-analysis script that performs both gene-level and genome-level tests. It should either be moved to the statistics directory as the active compute_statistics.R implementation or replaced by the separate gene-feature and genome-statistics scripts.

suppressPackageStartupMessages({
    library(argparse)
    library(dplyr)
    library(readr)
    library(tidyr)
    library(purrr)
    library(lme4)
    library(broom.mixed)
    library(tibble)
})

parser <- ArgumentParser(description = "Compute reproducible tAIpipe statistics")
parser$add_argument("--gene-table", required = TRUE)
parser$add_argument("--genome-table", required = TRUE)
parser$add_argument("--gene-feature-output", required = TRUE)
parser$add_argument("--genome-group-output", required = TRUE)
parser$add_argument(
    "--binary-features",
    default = "signal_peptide_present,tm_present,lcr_present,pfam_present"
)
parser$add_argument(
    "--genome-metrics",
    default = "mean_tAI,median_tAI,mean_GC3s,median_delta_ENC"
)
parser$add_argument("--group-variables", default = "phylum,lifestyle")
parser$add_argument("--fdr-method", default = "BH")
args <- parser$parse_args()

split_argument <- function(x) {
    values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    values[nzchar(values)]
}

safe_zscore <- function(x) {
    spread <- sd(x, na.rm = TRUE)
    if (!is.finite(spread) || spread == 0) {
        return(rep(NA_real_, length(x)))
    }
    as.numeric((x - mean(x, na.rm = TRUE)) / spread)
}

as_binary <- function(x) {
    if (is.logical(x)) {
        return(as.integer(x))
    }
    if (is.numeric(x)) {
        valid <- unique(x[!is.na(x)])
        if (all(valid %in% c(0, 1))) {
            return(as.integer(x))
        }
    }

    normalized <- tolower(trimws(as.character(x)))
    result <- rep(NA_integer_, length(normalized))
    result[normalized %in% c("1", "true", "t", "yes", "y", "present")] <- 1L
    result[normalized %in% c("0", "false", "f", "no", "n", "absent")] <- 0L
    result
}

fit_binary_feature <- function(gene_data, feature_name) {
    if (!feature_name %in% names(gene_data)) {
        return(tibble(
            analysis = "gene_level_mixed_model",
            feature = feature_name,
            term = NA_character_,
            estimate = NA_real_,
            std_error = NA_real_,
            statistic = NA_real_,
            p_value = NA_real_,
            conf_low = NA_real_,
            conf_high = NA_real_,
            n_genes = 0L,
            n_genomes = 0L,
            status = "missing_column"
        ))
    }

    model_data <- gene_data %>%
        transmute(
            sample,
            tAI_z,
            feature_value = as_binary(.data[[feature_name]])
        ) %>%
        filter(!is.na(sample), is.finite(tAI_z), !is.na(feature_value))

    n_genomes <- n_distinct(model_data$sample)
    if (
        nrow(model_data) < 20L ||
            n_genomes < 3L ||
            n_distinct(model_data$feature_value) != 2L
    ) {
        return(tibble(
            analysis = "gene_level_mixed_model",
            feature = feature_name,
            term = "feature_value",
            estimate = NA_real_,
            std_error = NA_real_,
            statistic = NA_real_,
            p_value = NA_real_,
            conf_low = NA_real_,
            conf_high = NA_real_,
            n_genes = nrow(model_data),
            n_genomes = n_genomes,
            status = "insufficient_data"
        ))
    }

    model <- tryCatch(
        lmer(tAI_z ~ feature_value + (1 | sample), data = model_data, REML = FALSE),
        error = function(error) error
    )

    if (inherits(model, "error")) {
        return(tibble(
            analysis = "gene_level_mixed_model",
            feature = feature_name,
            term = "feature_value",
            estimate = NA_real_,
            std_error = NA_real_,
            statistic = NA_real_,
            p_value = NA_real_,
            conf_low = NA_real_,
            conf_high = NA_real_,
            n_genes = nrow(model_data),
            n_genomes = n_genomes,
            status = paste0("model_error: ", conditionMessage(model))
        ))
    }

    tidy(model, effects = "fixed", conf.int = TRUE) %>%
        filter(term == "feature_value") %>%
        transmute(
            analysis = "gene_level_mixed_model",
            feature = feature_name,
            term,
            estimate,
            std_error = std.error,
            statistic,
            p_value = p.value,
            conf_low = conf.low,
            conf_high = conf.high,
            n_genes = nrow(model_data),
            n_genomes = n_genomes,
            status = ifelse(isSingular(model), "singular_fit", "ok")
        )
}

run_genome_group_test <- function(genome_data, metric_name, group_name) {
    if (!all(c(metric_name, group_name) %in% names(genome_data))) {
        return(tibble(
            analysis = "genome_level_group_test",
            metric = metric_name,
            group_variable = group_name,
            test = NA_character_,
            comparison = NA_character_,
            effect = NA_real_,
            statistic = NA_real_,
            p_value = NA_real_,
            n_genomes = 0L,
            n_groups = 0L,
            status = "missing_column"
        ))
    }

    test_data <- genome_data %>%
        transmute(
            value = as.numeric(.data[[metric_name]]),
            group = as.factor(.data[[group_name]])
        ) %>%
        filter(is.finite(value), !is.na(group)) %>%
        droplevels()

    group_count <- nlevels(test_data$group)
    if (nrow(test_data) < 4L || group_count < 2L) {
        return(tibble(
            analysis = "genome_level_group_test",
            metric = metric_name,
            group_variable = group_name,
            test = NA_character_,
            comparison = NA_character_,
            effect = NA_real_,
            statistic = NA_real_,
            p_value = NA_real_,
            n_genomes = nrow(test_data),
            n_groups = group_count,
            status = "insufficient_data"
        ))
    }

    if (group_count == 2L) {
        levels_group <- levels(test_data$group)
        result <- wilcox.test(value ~ group, data = test_data, exact = FALSE)
        medians <- test_data %>%
            group_by(group) %>%
            summarise(median = median(value), .groups = "drop")
        effect <- medians$median[[2]] - medians$median[[1]]

        return(tibble(
            analysis = "genome_level_group_test",
            metric = metric_name,
            group_variable = group_name,
            test = "Wilcoxon rank-sum",
            comparison = paste(levels_group, collapse = " vs "),
            effect = effect,
            statistic = unname(result$statistic),
            p_value = result$p.value,
            n_genomes = nrow(test_data),
            n_groups = group_count,
            status = "ok"
        ))
    }

    result <- kruskal.test(value ~ group, data = test_data)
    h_statistic <- unname(result$statistic)
    epsilon_squared <- max(
        0,
        (h_statistic - group_count + 1) / (nrow(test_data) - group_count)
    )

    tibble(
        analysis = "genome_level_group_test",
        metric = metric_name,
        group_variable = group_name,
        test = "Kruskal-Wallis",
        comparison = "global",
        effect = epsilon_squared,
        statistic = h_statistic,
        p_value = result$p.value,
        n_genomes = nrow(test_data),
        n_groups = group_count,
        status = "ok"
    )
}

gene_table <- read_tsv(args$gene_table, show_col_types = FALSE)
genome_table <- read_tsv(args$genome_table, show_col_types = FALSE)

required_gene_columns <- c("sample", "tAI")
missing_gene_columns <- setdiff(required_gene_columns, names(gene_table))
if (length(missing_gene_columns) > 0L) {
    stop("Gene table lacks columns: ", paste(missing_gene_columns, collapse = ", "))
}

gene_table <- gene_table %>%
    group_by(sample) %>%
    mutate(tAI_z = safe_zscore(tAI)) %>%
    ungroup()

binary_features <- split_argument(args$binary_features)
gene_results <- map_dfr(binary_features, ~ fit_binary_feature(gene_table, .x)) %>%
    mutate(q_value = p.adjust(p_value, method = args$fdr_method))

genome_metrics <- split_argument(args$genome_metrics)
group_variables <- split_argument(args$group_variables)
genome_combinations <- crossing(
    metric = genome_metrics,
    group_variable = group_variables
)

genome_results <- pmap_dfr(
    genome_combinations,
    ~ run_genome_group_test(genome_table, ..1, ..2)
) %>%
    mutate(q_value = p.adjust(p_value, method = args$fdr_method))

dir.create(dirname(args$gene_feature_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$genome_group_output), recursive = TRUE, showWarnings = FALSE)
write_tsv(gene_results, args$gene_feature_output, na = "NA")
write_tsv(genome_results, args$genome_group_output, na = "NA")
