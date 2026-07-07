"""
Tests whether Gene Ontology terms are enriched or depleted among genes with high or low tAI. It runs a genome-stratified Cochran–Mantel–Haenszel analysis and writes a result table containing common odds ratios, confidence intervals, genome counts, raw p-values, and adjusted p-values.
"""

rule compute_go_enrichment_cmh:
    input:
        gene_features=config["paths"]["gene_features"],
        go_dictionary=config["go_dictionary"]["path"],

    output:
        table=config["paths"]["go_enrichment"],

    params:
        tail_fraction=config["go_enrichment"]["tail_fraction"],
        max_tail_genes=config["go_enrichment"]["max_tail_genes"],
        min_genomes_with_term=config["go_enrichment"]["min_genomes_with_term"],
        min_informative_genomes=config["go_enrichment"]["min_informative_genomes"],
        min_total_genes_with_term=config["go_enrichment"]["min_total_genes_with_term"],
        fdr_method=config["go_enrichment"]["fdr_method"],

    log:
        f"{LOGS}/statistics/compute_go_enrichment_cmh.log"

    conda:
        "../envs/r.yaml"

    message:
        "Computing GO enrichment with genome-stratified CMH tests"

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
            --fdr-method {params.fdr_method:q} \
            > {log:q} 2>&1
        """
