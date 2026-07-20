"""
Runs the statistical analyses defined for the canonical gene and genome tables. It passes the configured features, metrics, grouping variables, and multiple-testing method to the analysis script and writes separate result tables for gene-level effects and genome-level group comparisons.
"""

rule compute_statistics:
    input:
        gene_table=config["paths"]["gene_features"],
        genome_table=config["paths"]["genome_summary"]
        
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
                ["mean_tAI", "mean_CAI", "mean_GC", "mean_GC3s", "mean_ENC", "mean_delta_ENC"],
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


rule compute_gene_feature_effects_per_genome:
    input:
        gene_table=config["paths"]["gene_features"]

    output:
        per_genome=f"{DATA_STATISTICS}/gene_feature_per_genome_effects.tsv",
        meta=f"{DATA_STATISTICS}/gene_feature_effect_meta_tests.tsv"

    params:
        features=",".join(config["statistics"]["binary_features"]),
        covariates=",".join(config["statistics"].get("gene_covariates", [])),
        min_genes=config.get("statistics", {}).get("per_genome_feature_effects", {}).get("min_genes", 100),
        min_class_genes=config.get("statistics", {}).get("per_genome_feature_effects", {}).get("min_class_genes", 10),
        min_genomes=config.get("statistics", {}).get("per_genome_feature_effects", {}).get("min_genomes", 5),
        fdr_method=config["statistics"].get("fdr_method", "BH")

    log:
        f"{LOGS}/statistics/compute_gene_feature_effects_per_genome.log"

    conda:
        "../envs/r_statistics.yaml"

    shell:
        """
        set -euo pipefail

        mkdir -p "$(dirname {output.per_genome:q})" "$(dirname {log:q})"

        Rscript workflow/scripts/statistics/compute_gene_feature_effects_per_genome.R \
          --gene-table {input.gene_table:q} \
          --features {params.features:q} \
          --covariates {params.covariates:q} \
          --min-genes {params.min_genes} \
          --min-class-genes {params.min_class_genes} \
          --min-genomes {params.min_genomes} \
          --fdr-method {params.fdr_method:q} \
          --per-genome-output {output.per_genome:q} \
          --meta-output {output.meta:q} \
          > {log:q} 2>&1
        """


rule compute_lifestyle_within_phylum_permutations:
    input:
        genome_table=config["paths"]["genome_summary"]

    output:
        f"{DATA_STATISTICS}/lifestyle_within_phylum_permutation_tests.tsv"

    params:
        responses=",".join(config.get("statistics", {}).get("lifestyle_within_phylum_permutations", {}).get(
            "responses",
            ["mean_tAI", "mean_CAI", "mean_GC", "mean_GC3s", "mean_ENC", "mean_delta_ENC"],
        )),
        predictor=config.get("statistics", {}).get("lifestyle_within_phylum_permutations", {}).get("predictor", "lifestyle"),
        strata=config.get("statistics", {}).get("lifestyle_within_phylum_permutations", {}).get("strata", "phylum"),
        n_perm=config.get("statistics", {}).get("lifestyle_within_phylum_permutations", {}).get("n_perm", 9999),
        seed=config.get("statistics", {}).get("lifestyle_within_phylum_permutations", {}).get("seed", 1),
        fdr_method=config.get("statistics", {}).get("fdr_method", "BH")

    log:
        f"{LOGS}/statistics/lifestyle_within_phylum_permutation_tests.log"

    conda:
        "../envs/r_statistics.yaml"

    shell:
        """
        set -euo pipefail

        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"

        Rscript workflow/scripts/statistics/restricted_lifestyle_permutations.R \
          --genome-table {input.genome_table:q} \
          --responses {params.responses:q} \
          --predictor {params.predictor:q} \
          --strata {params.strata:q} \
          --n-perm {params.n_perm} \
          --seed {params.seed} \
          --fdr-method {params.fdr_method:q} \
          --output {output:q} \
          > {log:q} 2>&1
        """