# Replace only rule plot_codon_profile_overview in workflow/rules/12_plots.smk
# with this version if your plotting script reads results/tables/codon_profiles.tsv
# directly.

rule plot_codon_profile_overview:
    input:
        codon_profiles=config["paths"]["codon_profiles"],
        plotting_config=PLOTTING_CONFIG,
        plotting_utils=R_PLOTTING_UTILS,
        table_validation_utils=R_TABLE_VALIDATION_UTILS,
        plot_data_utils=R_PLOT_DATA_UTILS,
        label_utils=R_LABEL_UTILS,

    output:
        variability_png=f"{DATA_PLOTS}/codon_profiles/codon_usage_variability.png",
        variability_pdf=f"{DATA_PLOTS}/codon_profiles/codon_usage_variability.pdf",
        method=f"{DATA_PLOTS}/codon_profiles/codon_usage_variability_method.tsv",
        heatmap_phylum=f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_phylum_large_codon_x.png",
        heatmap_lifestyle=f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_lifestyle_large_codon_x.png",

    params:
        output_dir=f"{DATA_PLOTS}/codon_profiles",
        formats=config.get("plots", {}).get("output_formats", ["png", "pdf"]),
        heatmap_value=config.get("codon_profile_plots", {}).get("heatmap_value", "auto"),
        usage_value=config.get("codon_profile_plots", {}).get("usage_value", "genome_RSCU"),
        top_n_variable_codons=config.get("codon_profile_plots", {}).get("top_n_variable_codons", 30),
        heatmap_large_width=config.get("codon_profile_plots", {}).get("heatmap_large_width", 16.5),
        heatmap_large_height=config.get("codon_profile_plots", {}).get("heatmap_large_height", 23.4),
        heatmap_min_group_n=config.get("codon_profile_plots", {}).get("heatmap_min_group_n", 5),

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

    shell:
        """
        set -euo pipefail

        mkdir -p {params.output_dir:q} "$(dirname {log.stdout:q})"

        Rscript workflow/scripts/plots/plot_codon_profile_overview.R \
            --codon-profiles {input.codon_profiles:q} \
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
