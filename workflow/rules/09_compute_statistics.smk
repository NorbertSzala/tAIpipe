"""
Runs the statistical analyses defined for the canonical gene and genome tables. It passes the configured features, metrics, grouping variables, and multiple-testing method to the analysis script and writes separate result tables for gene-level effects and genome-level group comparisons.
"""

rule compute_statistics:
    input:
        gene_table=f"{DATA_TABLES}/gene_features.tsv",
        genome_table=f"{DATA_TABLES}/genome_summary.tsv",

    output:
        gene_tests = config['paths']['gene_feature_tests'],
        genome_tests = config['paths']['genome_group_tests']

    params:
        binary_features=",".join(
            config.get("statistics", {}).get(
                "binary_features",
                [
                    "signal_peptide_present",
                    "tm_present",
                    "lcr_present",
                    "pfam_present",
                ],
            )
        ),
        gene_covariates = ",".join(
            config.get('statistics', {}).get('gene_covariates', ['log_protein_length_aa', "GC3s"],)
        ),

        genome_metrics=",".join(
            config.get("statistics", {}).get(
                "genome_metrics",
                ["mean_tAI", "median_tAI", "mean_GC3s", "median_delta_ENC"],
            )
        ),

        group_variables=",".join(
            config.get("statistics", {}).get(
                "group_variables",
                ["phylum", "lifestyle"],
            )
        ),
        fdr_method=config.get("statistics", {}).get("fdr_method", "BH"),

    log:
        f"{LOGS}/aggregated/compute_statistics.log"

    message:
        "Counting genome statistics"

    conda:
        "../envs/r.yaml"

    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.gene_tests})" "$(dirname {log})"

        Rscript workflow/scripts/statistics/compute_statistics.R \
          --gene-table {input.gene_table} \
          --genome-table {input.genome_table} \
          --gene-feature-output {output.gene_tests} \
          --genome-group-output {output.genome_tests} \
          --binary-features {params.binary_features} \
          --gene-covariates {params.gene_covariates} \
          --genome-metrics {params.genome_metrics} \
          --group-variables {params.group_variables} \
          --fdr-method {params.fdr_method} \
          > {log} 2>&1
        """
