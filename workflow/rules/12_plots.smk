# ------------------------------------------------------------------------------
# --- Static plotting rules -----------------------------------------------------
# ------------------------------------------------------------------------------
# These rules generate publication/report-ready static figures from plot-ready
# summary tables. Plotting scripts should not recompute core statistics.
# ------------------------------------------------------------------------------

PLOTTING_CONFIG = "config/plotting.yaml"

R_PLOTTING_ENV = "../envs/r_plotting.yaml"

CODON_HEATMAP_MIN_GROUP_N = config.get("codon_profile_plots", {}).get(
    "heatmap_min_group_n",
    5,
)

R_PLOTTING_UTILS = workflow.source_path("../scripts/lib/plotting_utils.R")
R_TABLE_VALIDATION_UTILS = workflow.source_path("../scripts/lib/table_validation_utils.R")
R_PLOT_DATA_UTILS = workflow.source_path("../scripts/lib/plot_data_utils.R")
R_LABEL_UTILS = workflow.source_path("../scripts/lib/label_utils.R")
R_PLOT_STYLE_HELPERS = workflow.source_path("../scripts/lib/plot_style_helpers.R")


# Plot tAI-tail feature summaries.
# Produces compact figures comparing all genes with upper/lower tAI tails.
rule plot_tai_extreme_features:
    input:
        binary_summary=rules.summarize_tai_extreme_features.output.binary_summary,
        continuous_summary=rules.summarize_tai_extreme_features.output.continuous_summary,
        plotting_config=PLOTTING_CONFIG,
        plotting_utils=R_PLOTTING_UTILS,
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        binary_png=f"{DATA_PLOTS}/tai_extremes/binary_features.png",
        binary_pdf=f"{DATA_PLOTS}/tai_extremes/binary_features.pdf",
        continuous_png=f"{DATA_PLOTS}/tai_extremes/continuous_features.png",
        continuous_pdf=f"{DATA_PLOTS}/tai_extremes/continuous_features.pdf"

    params:
        plot_group="tai_extremes",
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"])

    log:
        stdout=f"{LOGS}/plots/plot_tai_extreme_features.stdout.log",
        stderr=f"{LOGS}/plots/plot_tai_extreme_features.stderr.log"

    benchmark:
        f"{BENCHMARKS}/plots/plot_tai_extreme_features.tsv"

    threads: 1

    resources:
        mem_mb=4000,
        runtime=30

    conda:
        R_PLOTTING_ENV

    message:
        "Rendering static plots for tAI-tail binary and continuous feature summaries."

    script:
        "../scripts/plots/plot_tai_extreme_features.R"


# Plot gene-level feature overview.
# Uses precomputed distribution/effect summaries and does not run statistical
# tests directly inside the plotting script.
rule plot_gene_feature_overview:
    input:
        distribution_summary=rules.summarize_gene_feature_plots.output.distribution_summary,
        effect_summary=rules.summarize_gene_feature_plots.output.effect_summary,
        plotting_config=PLOTTING_CONFIG,
        plotting_utils=R_PLOTTING_UTILS,
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        binary_png=f"{DATA_PLOTS}/gene_features/tai_by_binary_features.png",
        binary_pdf=f"{DATA_PLOTS}/gene_features/tai_by_binary_features.pdf",
        continuous_png=f"{DATA_PLOTS}/gene_features/tai_vs_continuous_features.png",
        continuous_pdf=f"{DATA_PLOTS}/gene_features/tai_vs_continuous_features.pdf",
        effects_png=f"{DATA_PLOTS}/gene_features/gene_feature_effects.png",
        effects_pdf=f"{DATA_PLOTS}/gene_features/gene_feature_effects.pdf"

    params:
        plot_group="gene_features",
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"]),
        max_categories=config.get("plots", {}).get("max_categories_without_faceting", 8)

    log:
        stdout=f"{LOGS}/plots/plot_gene_feature_overview.stdout.log",
        stderr=f"{LOGS}/plots/plot_gene_feature_overview.stderr.log"

    benchmark:
        f"{BENCHMARKS}/plots/plot_gene_feature_overview.tsv"

    threads: 1

    resources:
        mem_mb=5000,
        runtime=45

    conda:
        R_PLOTTING_ENV

    message:
        "Rendering static gene-level feature overview plots."

    script:
        "../scripts/plots/plot_gene_feature_overview.R"


# Plot genome-level metric overview.
# This rule visualizes genome-level summaries by metadata groups such as phylum
# and lifestyle.
rule plot_genome_metric_overview:
    input:
        metric_summary=rules.summarize_genome_metric_plots.output.metric_summary,
        effect_summary=rules.summarize_genome_metric_plots.output.effect_summary,
        plotting_config=PLOTTING_CONFIG,
        plotting_utils=R_PLOTTING_UTILS,
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        phylum_png=f"{DATA_PLOTS}/genome_metrics/genome_metrics_by_phylum.png",
        phylum_pdf=f"{DATA_PLOTS}/genome_metrics/genome_metrics_by_phylum.pdf",
        lifestyle_png=f"{DATA_PLOTS}/genome_metrics/genome_metrics_by_lifestyle.png",
        lifestyle_pdf=f"{DATA_PLOTS}/genome_metrics/genome_metrics_by_lifestyle.pdf"

    params:
        plot_group="genome_metrics",
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"]),
        top_n_categories=config.get("plots", {}).get("default_top_n_categories", 12)

    log:
        stdout=f"{LOGS}/plots/plot_genome_metric_overview.stdout.log",
        stderr=f"{LOGS}/plots/plot_genome_metric_overview.stderr.log"

    benchmark:
        f"{BENCHMARKS}/plots/plot_genome_metric_overview.tsv"

    threads: 1

    resources:
        mem_mb=5000,
        runtime=45

    conda:
        R_PLOTTING_ENV

    message:
        "Rendering static genome-level metric plots by phylum and lifestyle."

    script:
        "../scripts/plots/plot_genome_metric_overview.R"


# Plot codon and tRNA profile overview.
# This rule reads the canonical codon_profiles.tsv directly because the current
# plotting script builds both its plot-specific summaries and figures in one pass.
rule plot_codon_profile_overview:
    input:
        codon_profiles=config["paths"]["codon_profiles"],
        genome_summary=config["paths"]["genome_summary"],
        plotting_config=PLOTTING_CONFIG,
        plotting_utils=R_PLOTTING_UTILS,
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        variability_png=f"{DATA_PLOTS}/codon_profiles/codon_usage_variability.png",
        variability_pdf=f"{DATA_PLOTS}/codon_profiles/codon_usage_variability.pdf",
        method=f"{DATA_PLOTS}/codon_profiles/codon_usage_variability_method.tsv",
        heatmap_phylum_png=f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_phylum_large_codon_x.png",
        heatmap_phylum_pdf=f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_phylum_large_codon_x.pdf",
        heatmap_lifestyle_png=f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_lifestyle_large_codon_x.png",
        heatmap_lifestyle_pdf=f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_lifestyle_large_codon_x.pdf",
        heatmap_phylum_collapsed_png=f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_phylum_count{CODON_HEATMAP_MIN_GROUP_N}_large_codon_x.png",
        heatmap_phylum_collapsed_pdf=f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_phylum_count{CODON_HEATMAP_MIN_GROUP_N}_large_codon_x.pdf",
        heatmap_lifestyle_collapsed_png=f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_lifestyle_count{CODON_HEATMAP_MIN_GROUP_N}_large_codon_x.png",
        heatmap_lifestyle_collapsed_pdf=f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_lifestyle_count{CODON_HEATMAP_MIN_GROUP_N}_large_codon_x.pdf",
        dendrogram_phylum_png=f"{DATA_PLOTS}/codon_profiles/trna_weights_dendrogram_heatmap_by_phylum.png",
        dendrogram_phylum_pdf=f"{DATA_PLOTS}/codon_profiles/trna_weights_dendrogram_heatmap_by_phylum.pdf",
        dendrogram_lifestyle_png=f"{DATA_PLOTS}/codon_profiles/trna_weights_dendrogram_heatmap_by_lifestyle.png",
        dendrogram_lifestyle_pdf=f"{DATA_PLOTS}/codon_profiles/trna_weights_dendrogram_heatmap_by_lifestyle.pdf",
        dendrogram_method=f"{DATA_PLOTS}/codon_profiles/trna_weights_dendrogram_method.tsv",
        annotation_association=f"{DATA_PLOTS}/codon_profiles/trna_weights_annotation_association.tsv"

    params:
        output_dir=f"{DATA_PLOTS}/codon_profiles",
        formats=",".join(config.get("plots", {}).get("output_formats", ["png", "pdf"])),
        heatmap_value=config.get("codon_profile_plots", {}).get("heatmap_value", "auto"),
        usage_value=config.get("codon_profile_plots", {}).get("usage_value", "genome_RSCU"),
        top_n_variable_codons=config.get("codon_profile_plots", {}).get("top_n_variable_codons", 30),
        heatmap_large_width=config.get("codon_profile_plots", {}).get("heatmap_large_width", 16.5),
        heatmap_large_height=config.get("codon_profile_plots", {}).get("heatmap_large_height", 23.4),
        heatmap_min_group_n=CODON_HEATMAP_MIN_GROUP_N

    log:
        stdout=f"{LOGS}/plots/plot_codon_profile_overview.stdout.log",
        stderr=f"{LOGS}/plots/plot_codon_profile_overview.stderr.log"

    benchmark:
        f"{BENCHMARKS}/plots/plot_codon_profile_overview.tsv"

    threads: 1

    resources:
        mem_mb=8000,
        runtime=60

    conda:
        R_PLOTTING_ENV

    message:
        "Rendering static codon/tRNA profile heatmaps and variability plots."

    shell:
        """
        set -euo pipefail

        mkdir -p {params.output_dir:q} "$(dirname {log.stdout:q})" "$(dirname {log.stderr:q})"

        Rscript workflow/scripts/plots/plot_codon_profile_overview.R \
            --codon-profiles {input.codon_profiles:q} \
            --genome-summary {input.genome_summary:q} \
            --output-dir {params.output_dir:q} \
            --formats {params.formats:q} \
            --heatmap-value {params.heatmap_value:q} \
            --usage-value {params.usage_value:q} \
            --top-n-variable-codons {params.top_n_variable_codons} \
            --heatmap-large-width {params.heatmap_large_width} \
            --heatmap-large-height {params.heatmap_large_height} \
            --heatmap-min-group-n {params.heatmap_min_group_n} \
            > {log.stdout:q} 2> {log.stderr:q}
        """


rule plot_go_enrichment_overview:
    input:
        top_terms=rules.summarize_go_enrichment_plots.output.top_terms,
        namespace_summary=rules.summarize_go_enrichment_plots.output.namespace_summary,
        plotting_config=PLOTTING_CONFIG,
        plotting_utils=R_PLOTTING_UTILS,
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        dotplot_png=f"{DATA_PLOTS}/go_enrichment/go_top_terms_dotplot.png",
        dotplot_pdf=f"{DATA_PLOTS}/go_enrichment/go_top_terms_dotplot.pdf",
        forest_png=f"{DATA_PLOTS}/go_enrichment/go_odds_ratio_forest.png",
        forest_pdf=f"{DATA_PLOTS}/go_enrichment/go_odds_ratio_forest.pdf"

    params:
        plot_group="go_enrichment",
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"]),
        top_n=config.get("plots", {}).get("default_top_n_terms", 20)

    log:
        stdout=f"{LOGS}/plots/plot_go_enrichment_overview.stdout.log",
        stderr=f"{LOGS}/plots/plot_go_enrichment_overview.stderr.log"

    benchmark:
        f"{BENCHMARKS}/plots/plot_go_enrichment_overview.tsv"

    threads: 1

    resources:
        mem_mb=5000,
        runtime=45

    conda:
        R_PLOTTING_ENV

    message:
        "Rendering static GO enrichment dotplot and forest plot."

    script:
        "../scripts/plots/plot_go_enrichment_overview.R"


# Plot gene-level and genome-level correlation summaries.
# The two matrices are kept separate because genes and genomes are different
# biological/statistical units.
rule plot_correlation_overview:
    input:
        gene_correlations=rules.summarize_correlation_plots.output.gene_correlations,
        genome_correlations=rules.summarize_correlation_plots.output.genome_correlations,
        plotting_config=PLOTTING_CONFIG,
        plotting_utils=R_PLOTTING_UTILS,
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        gene_png=f"{DATA_PLOTS}/correlations/gene_level_correlation_matrix.png",
        gene_pdf=f"{DATA_PLOTS}/correlations/gene_level_correlation_matrix.pdf",
        genome_png=f"{DATA_PLOTS}/correlations/genome_level_correlation_matrix.png",
        genome_pdf=f"{DATA_PLOTS}/correlations/genome_level_correlation_matrix.pdf"

    params:
        plot_group="correlations",
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"])

    log:
        stdout=f"{LOGS}/plots/plot_correlation_overview.stdout.log",
        stderr=f"{LOGS}/plots/plot_correlation_overview.stderr.log"

    benchmark:
        f"{BENCHMARKS}/plots/plot_correlation_overview.tsv"

    threads: 1

    resources:
        mem_mb=5000,
        runtime=45

    conda:
        R_PLOTTING_ENV

    message:
        "Rendering static gene-level and genome-level correlation heatmaps."

    script:
        "../scripts/plots/plot_correlation_overview.R"


# Plot technical QC overview.
# This rule creates diagnostic figures from per-genome tRNA/tAI QC summaries.
rule plot_qc_overview:
    input:
        qc_summary=rules.summarize_qc_plots.output.qc_summary,
        plotting_config=PLOTTING_CONFIG,
        plotting_utils=R_PLOTTING_UTILS,
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS

    output:
        qc_png=f"{DATA_PLOTS}/qc/trna_qc_flags.png",
        qc_pdf=f"{DATA_PLOTS}/qc/trna_qc_flags.pdf"

    params:
        plot_group="qc",
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"])

    log:
        stdout=f"{LOGS}/plots/plot_qc_overview.stdout.log",
        stderr=f"{LOGS}/plots/plot_qc_overview.stderr.log"

    benchmark:
        f"{BENCHMARKS}/plots/plot_qc_overview.tsv"

    threads: 1

    resources:
        mem_mb=4000,
        runtime=30

    conda:
        R_PLOTTING_ENV

    message:
        "Rendering static tRNA/tAI QC overview plots."

    script:
        "../scripts/plots/plot_qc_overview.R"
