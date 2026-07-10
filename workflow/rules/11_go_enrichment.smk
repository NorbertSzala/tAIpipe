"""
GO enrichment rules.

Computes enrichment/depletion of GO terms in genes from high- and low-tAI tails
using a genome-stratified Cochran-Mantel-Haenszel analysis. The output path is
controlled by config["paths"]["go_enrichment"].
"""

rule compute_go_enrichment_cmh:
    input:
        gene_features=config["paths"]["gene_features"],
        go_dictionary=config["go_dictionary"]["path"],

    output:
        table=config["paths"]["go_enrichment"],

    params:
        tail_fraction=config.get("go_enrichment", {}).get("tail_fraction", 0.10),
        max_tail_genes=config.get("go_enrichment", {}).get("max_tail_genes", 0),
        min_genomes_with_term=config.get("go_enrichment", {}).get("min_genomes_with_term", 5),
        min_informative_genomes=config.get("go_enrichment", {}).get("min_informative_genomes", 3),
        min_total_genes_with_term=config.get("go_enrichment", {}).get("min_total_genes_with_term", 10),
        min_genes_with_go_total=config.get("go_enrichment", {}).get("min_genes_with_go_total", 100),
        min_genes_with_go_per_sample=config.get("go_enrichment", {}).get("min_genes_with_go_per_sample", 10),
        fdr_method=config.get("go_enrichment", {}).get("fdr_method", "BH"),

    log:
        f"{LOGS}/statistics/compute_go_enrichment_cmh.log"

    conda:
        "../envs/r.yaml"

    message:
        "Computing GO enrichment with genome-stratified CMH tests."

    shell:
        """
        set -euo pipefail

        mkdir -p \
            "$(dirname {output.table:q})" \
            "$(dirname {log:q})"

        Rscript workflow/scripts/statistics/compute_go_enrichment_cmh.R \
            --gene-features {input.gene_features:q} \
            --go-dictionary {input.go_dictionary:q} \
            --output {output.table:q} \
            --tail-fraction {params.tail_fraction} \
            --max-tail-genes {params.max_tail_genes} \
            --min-genomes-with-term {params.min_genomes_with_term} \
            --min-informative-genomes {params.min_informative_genomes} \
            --min-total-genes-with-term {params.min_total_genes_with_term} \
            --min-genes-with-go-total {params.min_genes_with_go_total} \
            --min-genes-with-go-per-sample {params.min_genes_with_go_per_sample} \
            --fdr-method {params.fdr_method:q} \
            > {log:q} 2>&1
        """
