# ------------------------------------------------------------------------------
# --- Plot-ready statistical summaries -----------------------------------------
# ------------------------------------------------------------------------------
# These rules convert canonical result tables and statistical outputs into
# compact, plot-ready TSV files. Plotting scripts should read these summaries
# instead of recomputing filters, tails, top-N terms, or aggregations.
# ------------------------------------------------------------------------------

PLOT_SUMMARIES = f"{DATA_STATISTICS}/plot_summaries"

R_STATISTICS_ENV = "../envs/r_statistics.yaml"

R_TABLE_VALIDATION_UTILS = workflow.source_path("../scripts/lib/table_validation_utils.R")
R_PLOT_DATA_UTILS = workflow.source_path("../scripts/lib/plot_data_utils.R")
R_LABEL_UTILS = workflow.source_path("../scripts/lib/label_utils.R")


# Summarize features among genes from the upper/lower tAI tails.
# Unit of interpretation: genes nested within genomes; tails are selected
# independently per sample to avoid domination by gene-rich genomes.
rule summarize_tai_extreme_features:
    input:
        gene_features=config["paths"]["gene_features"],
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        membership=f"{DATA_STATISTICS}/tai_extremes/tail_membership.tsv",
        binary_summary=f"{DATA_STATISTICS}/tai_extremes/binary_feature_summary.tsv",
        continuous_summary=f"{DATA_STATISTICS}/tai_extremes/continuous_feature_summary.tsv"

    params:
        tai_extremes=config["statistics"]["tai_extremes"]

    log:
        stdout=f"{LOGS}/statistics/summarize_tai_extreme_features.stdout.log",
        stderr=f"{LOGS}/statistics/summarize_tai_extreme_features.stderr.log"

    benchmark:
        f"{BENCHMARKS}/statistics/summarize_tai_extreme_features.tsv"

    threads: 1

    resources:
        mem_mb=4000,
        runtime=30

    conda:
        R_STATISTICS_ENV

    message:
        "Summarizing binary and continuous gene features in top/bottom tAI tails."

    script:
        "../scripts/statistics/summarize_tai_extreme_features.R"


# Prepare plot-ready summaries for gene-level feature associations.
# This rule should not perform the primary statistical tests; it reshapes and
# summarizes gene_features.tsv and gene_feature_tests.tsv for visualization.
rule summarize_gene_feature_plots:
    input:
        gene_features=config["paths"]["gene_features"],
        gene_tests=config["paths"]["gene_feature_tests"],
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        distribution_summary=f"{PLOT_SUMMARIES}/gene_feature_distribution_summary.tsv",
        effect_summary=f"{PLOT_SUMMARIES}/gene_feature_effect_summary.tsv"

    params:
        statistics=config["statistics"],
        max_categories=config.get("plots", {}).get("max_categories_without_faceting", 8)

    log:
        stdout=f"{LOGS}/statistics/summarize_gene_feature_plots.stdout.log",
        stderr=f"{LOGS}/statistics/summarize_gene_feature_plots.stderr.log"

    benchmark:
        f"{BENCHMARKS}/statistics/summarize_gene_feature_plots.tsv"

    threads: 1

    resources:
        mem_mb=4000,
        runtime=30

    conda:
        R_STATISTICS_ENV

    message:
        "Preparing plot-ready summaries for gene-level feature associations."

    script:
        "../scripts/statistics/summarize_gene_feature_plots.R"


# Prepare genome-level metric summaries for phylum/lifestyle comparisons.
# Unit of interpretation: genome. This keeps genome-level plots separate from
# gene-level plots and avoids pseudo-replication.
rule summarize_genome_metric_plots:
    input:
        genome_summary=config["paths"]["genome_summary"],
        genome_tests=config["paths"]["genome_group_tests"],
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        metric_summary=f"{PLOT_SUMMARIES}/genome_metric_summary.tsv",
        effect_summary=f"{PLOT_SUMMARIES}/genome_group_effect_summary.tsv"

    params:
        statistics=config["statistics"],
        top_n_categories=config.get("plots", {}).get("default_top_n_categories", 12)

    log:
        stdout=f"{LOGS}/statistics/summarize_genome_metric_plots.stdout.log",
        stderr=f"{LOGS}/statistics/summarize_genome_metric_plots.stderr.log"

    benchmark:
        f"{BENCHMARKS}/statistics/summarize_genome_metric_plots.tsv"

    threads: 1

    resources:
        mem_mb=4000,
        runtime=30

    conda:
        R_STATISTICS_ENV

    message:
        "Preparing plot-ready summaries for genome-level metrics and group effects."

    script:
        "../scripts/statistics/summarize_genome_metric_plots.R"


# Summarize codon/tRNA profiles for heatmaps and variability plots.
# This rule should reduce codon_profiles.tsv into compact matrices and
# long-format summaries suitable for static figures.
rule summarize_codon_profile_plots:
    input:
        codon_profiles=config["paths"]["codon_profiles"],
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        heatmap=f"{PLOT_SUMMARIES}/codon_profile_heatmap.tsv",
        variability=f"{PLOT_SUMMARIES}/codon_profile_variability.tsv",
        reference_comparison=f"{PLOT_SUMMARIES}/rscu_reference_comparison.tsv"

    params:
        top_n_categories=config.get("plots", {}).get("default_top_n_categories", 12)

    log:
        stdout=f"{LOGS}/statistics/summarize_codon_profile_plots.stdout.log",
        stderr=f"{LOGS}/statistics/summarize_codon_profile_plots.stderr.log"

    benchmark:
        f"{BENCHMARKS}/statistics/summarize_codon_profile_plots.tsv"

    threads: 1

    resources:
        mem_mb=6000,
        runtime=45

    conda:
        R_STATISTICS_ENV

    message:
        "Preparing codon-profile summaries for heatmaps and variability plots."

    script:
        "../scripts/statistics/summarize_codon_profile_plots.R"


# Select and format GO enrichment results for static plots.
# This rule should only work on the final enrichment table; enrichment statistics
# themselves should remain in the dedicated GO enrichment rule.
rule summarize_go_enrichment_plots:
    input:
        go_enrichment=config["paths"]["go_enrichment"],
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        top_terms=f"{PLOT_SUMMARIES}/go_enrichment_top_terms.tsv",
        namespace_summary=f"{PLOT_SUMMARIES}/go_enrichment_by_namespace.tsv"

    params:
        top_n=config.get("plots", {}).get("default_top_n_terms", 20),
        fdr_method=config.get("go_enrichment", {}).get("fdr_method", "BH")

    log:
        stdout=f"{LOGS}/statistics/summarize_go_enrichment_plots.stdout.log",
        stderr=f"{LOGS}/statistics/summarize_go_enrichment_plots.stderr.log"

    benchmark:
        f"{BENCHMARKS}/statistics/summarize_go_enrichment_plots.tsv"

    threads: 1

    resources:
        mem_mb=4000,
        runtime=30

    conda:
        R_STATISTICS_ENV

    message:
        "Preparing top GO enrichment terms and namespace summaries for plotting."

    script:
        "../scripts/statistics/summarize_go_enrichment_plots.R"


# Compute compact correlation tables used by correlation heatmaps.
# Keep gene-level and genome-level correlations separate because their biological
# units and interpretation are different.
rule summarize_correlation_plots:
    input:
        gene_features=config["paths"]["gene_features"],
        genome_summary=config["paths"]["genome_summary"],
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        gene_correlations=f"{PLOT_SUMMARIES}/gene_level_correlations.tsv",
        genome_correlations=f"{PLOT_SUMMARIES}/genome_level_correlations.tsv"

    params:
        gene_variables=config.get("plots", {}).get(
            "gene_correlation_variables",
            ["tAI", "CAI", "GC3s", "ENC", "protein_length_aa"]
        ),
        genome_variables=config.get("plots", {}).get(
            "genome_correlation_variables",
            ["mean_tAI", "median_tAI", "mean_GC3s", "median_delta_ENC"]
        ),
        method=config.get("plots", {}).get("correlation_method", "spearman")

    log:
        stdout=f"{LOGS}/statistics/summarize_correlation_plots.stdout.log",
        stderr=f"{LOGS}/statistics/summarize_correlation_plots.stderr.log"

    benchmark:
        f"{BENCHMARKS}/statistics/summarize_correlation_plots.tsv"

    threads: 1

    resources:
        mem_mb=6000,
        runtime=45

    conda:
        R_STATISTICS_ENV

    message:
        "Computing gene-level and genome-level correlation summaries for plotting."

    script:
        "../scripts/statistics/summarize_correlation_plots.R"


# Collect per-genome tRNA/tAI QC files into a single plot-ready table.
# This is a technical diagnostic layer, not a biological result layer.
rule summarize_qc_plots:
    input:
        qc_files=expand(
            f"{PER_GENOME}/{{sample}}/qc/{{sample}}_trna_profile_qc.tsv",
            sample=SAMPLES,
        ),
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        qc_summary=f"{PLOT_SUMMARIES}/trna_qc_summary.tsv"

    params:
        samples=SAMPLES

    log:
        stdout=f"{LOGS}/statistics/summarize_qc_plots.stdout.log",
        stderr=f"{LOGS}/statistics/summarize_qc_plots.stderr.log"

    benchmark:
        f"{BENCHMARKS}/statistics/summarize_qc_plots.tsv"

    threads: 1

    resources:
        mem_mb=3000,
        runtime=20

    conda:
        R_STATISTICS_ENV

    message:
        "Collecting per-genome tRNA profile QC summaries for plotting."

    script:
        "../scripts/statistics/summarize_qc_plots.R"