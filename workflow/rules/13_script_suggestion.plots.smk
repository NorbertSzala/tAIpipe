# ------------------------------------------------------------------------------
# Optional script_suggestions-compatible plot layer
# ------------------------------------------------------------------------------
# These rules recreate the biological meaning of historical exploratory plots
# from scripts_suggestions using canonical pipeline tables.
#
# They intentionally do not recompute tAI/CAI/ENC or enrichment statistics.
# Missing or ambiguous historical inputs are recorded as XXXXX_* files inside
# the corresponding plot directory.
# ------------------------------------------------------------------------------

SCRIPT_SUGGESTION_PLOT_ROOT = f"{DATA_PLOTS}/script_suggestions"


def script_suggestion_chosen_go_terms(wildcards):
    configured = config.get("script_suggestion_plots", {}).get(
        "chosen_go_terms",
        f"{DATA_STATISTICS}/chosen_GOterms.tsv",
    )
    if configured is None or str(configured).strip() in ["", "XXXXX", "auto"]:
        return f"{DATA_STATISTICS}/chosen_GOterms.tsv"
    return configured

rule plot_script_suggestion_feature_tai:
    input:
        gene_features=config["paths"]["gene_features"]
    output:
        outdir=directory(f"{SCRIPT_SUGGESTION_PLOT_ROOT}/features_tai")
    params:
        output_dir=f"{SCRIPT_SUGGESTION_PLOT_ROOT}/features_tai",
        metric=config.get("script_suggestion_plots", {}).get("feature_metric", "tAI"),
        top_fraction=config.get("script_suggestion_plots", {}).get("top_fraction", 0.10),
        min_genes_per_sample=config.get("statistics", {}).get("tai_extremes", {}).get("min_genes_per_sample", 100),
        min_tail_size=config.get("statistics", {}).get("tai_extremes", {}).get("min_tail_size", 20),
        max_violin_rows_per_group=config.get("script_suggestion_plots", {}).get("max_violin_rows_per_group", 50000),
        max_scatter_rows=config.get("script_suggestion_plots", {}).get("max_scatter_rows", 100000),
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"])
    conda:
        R_PLOTTING_ENV
    script:
        "../scripts/plots/plot_feature_tai_suggestion_set.R"


rule plot_script_suggestion_go_term_features:
    input:
        gene_features=config["paths"]["gene_features"],
        chosen_go_terms=script_suggestion_chosen_go_terms
    output:
        outdir=directory(f"{SCRIPT_SUGGESTION_PLOT_ROOT}/go_term_features")
    params:
        output_dir=f"{SCRIPT_SUGGESTION_PLOT_ROOT}/go_term_features",
        metric=config.get("script_suggestion_plots", {}).get("feature_metric", "tAI"),
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"])
    conda:
        R_PLOTTING_ENV
    script:
        "../scripts/plots/plot_go_term_feature_suggestion_set.R"


rule plot_script_suggestion_codon_trna:
    input:
        codon_profiles=config["paths"]["codon_profiles"]
    output:
        outdir=directory(f"{SCRIPT_SUGGESTION_PLOT_ROOT}/codon_trna")
    params:
        output_dir=f"{SCRIPT_SUGGESTION_PLOT_ROOT}/codon_trna",
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"])
    conda:
        R_PLOTTING_ENV
    script:
        "../scripts/plots/plot_codon_trna_suggestion_set.R"


rule plot_script_suggestion_genome_gc_tai:
    input:
        genome_summary=config["paths"]["genome_summary"]
    output:
        outdir=directory(f"{SCRIPT_SUGGESTION_PLOT_ROOT}/genome_gc_tai")
    params:
        output_dir=f"{SCRIPT_SUGGESTION_PLOT_ROOT}/genome_gc_tai",
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"])
    conda:
        R_PLOTTING_ENV
    script:
        "../scripts/plots/plot_genome_gc_tai_suggestion_set.R"


rule plot_script_suggestion_correlations:
    input:
        gene_features=config["paths"]["gene_features"],
        genome_summary=config["paths"]["genome_summary"]
    output:
        outdir=directory(f"{SCRIPT_SUGGESTION_PLOT_ROOT}/correlations")
    params:
        output_dir=f"{SCRIPT_SUGGESTION_PLOT_ROOT}/correlations",
        gene_variables=config.get("script_suggestion_plots", {}).get(
            "gene_correlation_variables",
            ["tAI", "CAI", "ENC", "GC", "GC3s", "protein_length_aa", "lcr_total_length", "tm_total_length"],
        ),
        genome_variables=config.get("script_suggestion_plots", {}).get(
            "genome_correlation_variables",
            ["mean_tAI", "median_tAI", "mean_CAI", "median_CAI", "mean_ENC", "median_ENC", "genome_gc", "mean_GC3s", "median_GC3s", "fraction_lcr", "fraction_tm", "fraction_signal_peptide", "fraction_pfam"],
        ),
        method=config.get("script_suggestion_plots", {}).get("correlation_method", "spearman"),
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"])
    conda:
        R_PLOTTING_ENV
    script:
        "../scripts/plots/plot_correlation_suggestion_set.R"
