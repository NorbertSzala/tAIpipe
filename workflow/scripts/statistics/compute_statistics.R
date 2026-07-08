#!/usr/bin/env Rscript
# Legacy combined statistical-analysis script that performs both gene-level and genome-level tests. It should either be moved to the statistics directory as the active compute_statistics.R implementation or replaced by the separate gene-feature and genome-statistics scripts.

suppressPackageStartupMessages({
    library(argparse)
    library(broom.mixed)
    library(dplyr)
    library(lme4)
    library(purrr)
    library(readr)
    library(tibble)
    library(tidyr)
})

parser <- ArgumentParser(description = "Compute tAIpipe statistical tests")
parser$add_argument("--gene-table", required = TRUE)
parser$add_argument("--genome-table", required = TRUE)
parser$add_argument("--gene-feature-output", required = TRUE)
parser$add_argument("--genome-group-output", required = TRUE)
parser$add_argument(
    "--binary-features",
    default = "signal_peptide_present,tm_present,lcr_present,pfam_present"
)
parser$add_argument(
    "--gene-covariates",
    default = "log_protein_length_aa,GC3s"
)
parser$add_argument(
    "--genome-metrics",
    default = "mean_tAI,median_tAI,mean_GC3s,median_delta_ENC"
)
parser$add_argument("--group-variables", default = "phylum,lifestyle")
parser$add_argument("--fdr-method", default = "BH")
args <- parser$parse_args()

split_argument <- function(x) {
    if (is.null(x) || is.na(x) || !nzchar(trimws(x))) {
        return(character())
    }
    values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    unique(values[nzchar(values)])
}

safe_zscore <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    finite <- is.finite(x)
    result <- rep(NA_real_, length(x))
    if (sum(finite) < 2L) {
        return(result)
    }
    spread <- sd(x[finite])
    if (!is.finite(spread) || spread == 0) {
        return(result)
    }
    result[finite] <- (x[finite] - mean(x[finite])) / spread
    result
}

as_binary <- function(x) {
    if (is.logical(x)) {
        return(as.integer(x))
    }
    if (is.numeric(x)) {
        valid <- unique(x[is.finite(x)])
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

empty_gene_result <- function(feature, status, n_genes = 0L, n_genomes = 0L,
                              covariates = character(), formula = NA_character_) {
    tibble(
        analysis = "gene_level_mixed_model",
        feature = feature,
        term = "feature_value",
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        n_genes = as.integer(n_genes),
        n_genomes = as.integer(n_genomes),
        covariates = paste(covariates, collapse = ";"),
        model_formula = formula,
        status = status
    )
}

fit_binary_feature <- function(gene_data, feature_name, requested_covariates) {
    if (!feature_name %in% names(gene_data)) {
        return(empty_gene_result(feature_name, "missing_feature_column"))
    }

    available_covariates <- intersect(requested_covariates, names(gene_data))
    model_data <- gene_data %>%
        transmute(
            sample = as.character(sample),
            tAI_z = suppressWarnings(as.numeric(tAI_z)),
            feature_value = as_binary(.data[[feature_name]]),
            across(all_of(available_covariates), ~ suppressWarnings(as.numeric(.x)))
        )

    # Covariates with no usable variation cannot be estimated and are omitted.
    usable_covariates <- available_covariates[vapply(
        model_data[available_covariates],
        function(x) {
            finite <- x[is.finite(x)]
            length(finite) >= 2L && length(unique(finite)) >= 2L
        },
        logical(1)
    )]

    keep_columns <- c("sample", "tAI_z", "feature_value", usable_covariates)
    model_data <- model_data %>%
        select(all_of(keep_columns)) %>%
        filter(
            !is.na(sample), nzchar(trimws(sample)),
            is.finite(tAI_z), !is.na(feature_value)
        )

    if (length(usable_covariates) > 0L) {
        complete_covariates <- complete.cases(model_data[usable_covariates]) &
            apply(model_data[usable_covariates], 1, function(row) all(is.finite(row)))
        model_data <- model_data[complete_covariates, , drop = FALSE]
    }

    n_genomes <- n_distinct(model_data$sample)
    fixed_terms <- c("feature_value", usable_covariates)
    formula_text <- paste(
        "tAI_z ~", paste(fixed_terms, collapse = " + "), "+ (1 | sample)"
    )

    if (nrow(model_data) < 20L || n_genomes < 3L ||
        n_distinct(model_data$feature_value) != 2L) {
        return(empty_gene_result(
            feature_name, "insufficient_data", nrow(model_data), n_genomes,
            usable_covariates, formula_text
        ))
    }

    model <- tryCatch(
        lmer(as.formula(formula_text), data = model_data, REML = FALSE),
        error = function(error) error
    )
    if (inherits(model, "error")) {
        return(empty_gene_result(
            feature_name,
            paste0("model_error: ", conditionMessage(model)),
            nrow(model_data), n_genomes, usable_covariates, formula_text
        ))
    }

    reduced_terms <- usable_covariates
    reduced_formula_text <- if (length(reduced_terms) > 0L) {
        paste("tAI_z ~", paste(reduced_terms, collapse = " + "), "+ (1 | sample)")
    } else {
        "tAI_z ~ 1 + (1 | sample)"
    }
    reduced_model <- tryCatch(
        lmer(as.formula(reduced_formula_text), data = model_data, REML = FALSE),
        error = function(error) error
    )
    if (inherits(reduced_model, "error")) {
        return(empty_gene_result(
            feature_name,
            paste0("reduced_model_error: ", conditionMessage(reduced_model)),
            nrow(model_data), n_genomes, usable_covariates, formula_text
        ))
    }

    likelihood_test <- tryCatch(
        anova(reduced_model, model, refit = FALSE),
        error = function(error) error
    )
    if (inherits(likelihood_test, "error")) {
        return(empty_gene_result(
            feature_name,
            paste0("likelihood_ratio_error: ", conditionMessage(likelihood_test)),
            nrow(model_data), n_genomes, usable_covariates, formula_text
        ))
    }
    p_value_lrt <- likelihood_test$`Pr(>Chisq)`[[2]]

    result <- tidy(
        model,
        effects = "fixed",
        conf.int = TRUE,
        conf.method = "Wald"
    ) %>%
        filter(term == "feature_value")
    if (nrow(result) != 1L) {
        return(empty_gene_result(
            feature_name, "feature_term_not_estimable", nrow(model_data), n_genomes,
            usable_covariates, formula_text
        ))
    }

    result %>%
        transmute(
            analysis = "gene_level_mixed_model",
            feature = feature_name,
            term,
            estimate,
            std_error = std.error,
            statistic,
            p_value = p_value_lrt,
            conf_low = conf.low,
            conf_high = conf.high,
            n_genes = nrow(model_data),
            n_genomes = n_genomes,
            covariates = paste(usable_covariates, collapse = ";"),
            model_formula = formula_text,
            status = ifelse(isSingular(model), "singular_fit", "ok")
        )
}

empty_genome_result <- function(metric, group, status, n = 0L, k = 0L) {
    tibble(
        analysis = "genome_level_group_test",
        metric = metric,
        group_variable = group,
        test = NA_character_,
        comparison = NA_character_,
        effect = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        n_genomes = as.integer(n),
        n_groups = as.integer(k),
        status = status
    )
}

run_genome_group_test <- function(genome_data, metric_name, group_name) {
    if (!all(c(metric_name, group_name) %in% names(genome_data))) {
        return(empty_genome_result(metric_name, group_name, "missing_column"))
    }

    test_data <- genome_data %>%
        transmute(
            value = suppressWarnings(as.numeric(.data[[metric_name]])),
            group = trimws(as.character(.data[[group_name]]))
        ) %>%
        filter(is.finite(value), !is.na(group), nzchar(group)) %>%
        mutate(group = factor(group)) %>%
        droplevels()

    group_count <- nlevels(test_data$group)
    if (nrow(test_data) < 4L || group_count < 2L) {
        return(empty_genome_result(
            metric_name, group_name, "insufficient_data", nrow(test_data), group_count
        ))
    }

    if (group_count == 2L) {
        if (any(table(test_data$group) < 2L)) {
            return(empty_genome_result(
                metric_name, group_name, "insufficient_group_size",
                nrow(test_data), group_count
            ))
        }
        levels_group <- levels(test_data$group)
        result <- wilcox.test(value ~ group, data = test_data, exact = FALSE)
        medians <- tapply(test_data$value, test_data$group, median)
        effect <- unname(medians[[2]] - medians[[1]])

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

gene_table <- read_tsv(args$gene_table, show_col_types = FALSE, progress = FALSE)
genome_table <- read_tsv(args$genome_table, show_col_types = FALSE, progress = FALSE)
required_gene <- c("sample", "tAI")
missing_gene <- setdiff(required_gene, names(gene_table))
if (length(missing_gene) > 0L) {
    stop("Gene table lacks columns: ", paste(missing_gene, collapse = ", "))
}

gene_table <- gene_table %>%
    group_by(sample) %>%
    mutate(tAI_z = safe_zscore(tAI)) %>%
    ungroup()

binary_features <- split_argument(args$binary_features)
gene_covariates <- split_argument(args$gene_covariates)
if (length(binary_features) == 0L) stop("No binary features were configured")

gene_results <- map_dfr(
    binary_features,
    ~ fit_binary_feature(gene_table, .x, gene_covariates)
) %>%
    mutate(q_value = p.adjust(p_value, method = args$fdr_method))

genome_metrics <- split_argument(args$genome_metrics)
group_variables <- split_argument(args$group_variables)
if (length(genome_metrics) == 0L || length(group_variables) == 0L) {
    stop("Genome metrics and group variables must not be empty")
}

genome_results <- crossing(
    metric = genome_metrics,
    group_variable = group_variables
) %>%
    pmap_dfr(~ run_genome_group_test(genome_table, ..1, ..2)) %>%
    mutate(q_value = p.adjust(p_value, method = args$fdr_method))

dir.create(dirname(args$gene_feature_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$genome_group_output), recursive = TRUE, showWarnings = FALSE)
write_tsv(gene_results, args$gene_feature_output, na = "NA")
write_tsv(genome_results, args$genome_group_output, na = "NA")
