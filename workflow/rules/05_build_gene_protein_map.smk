"""
Build protein-to-gene-to-CDS identifier mappings for every analysed genome.
"""

rule build_gene_protein_map:
    input:
        cds=get_cds,
        proteins=get_proteome,

    output:
        mapping=(
            f"{PER_GENOME}/{{sample}}/tables/"
            "gene_protein_map.tsv"
        ),

    params:
        sample=lambda wildcards: wildcards.sample,
        strict_proteome_match=False,
        max_missing_proteome_fraction=  config.get("gene_protein_map",{}).get("max_missing_proteome_fraction", 0.05),
        max_missing_proteome_count = config.get("gene_protein_map", {}).get("max_missing_proteome_count", 200)

    log:
        f"{LOGS}/{{sample}}/build_gene_protein_map.log"

    conda:
        "../envs/kofamscan.yaml"

    script:
        "../scripts/kofamscan/build_gene_protein_map.py"