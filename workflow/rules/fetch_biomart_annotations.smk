"""
Retrieves optional gene annotations from Ensembl BioMart for samples with configured BioMart identifiers. The rule produces a standardized annotation table that can be merged with the canonical gene-level dataset.

Requires explicit BioMart columns in metadata_dataset.tsv:
biomart_host, biomart_mart, biomart_dataset, biomart_id_filter,
biomart_id_attribute. """

def metadata_value(sample, column):
    value = samples_df.loc[sample, column]
    if pd.isna(value) or str(value).strip() == "":
        raise ValueError(f"Missing {column} for sample {sample}")
    return str(value).strip()

rule fetch_biomart_annotations:
    input:
        gene_table=BASE_GENE_TABLE,

    output:
        annotations=f"{PER_GENOME}/{{sample}}/annotations/{{sample}}_biomart.tsv",
    
    params:
        sample=lambda wildcards: wildcards.sample,
        host=lambda wildcards: metadata_value(wildcards.sample, "biomart_host"),
        mart=lambda wildcards: metadata_value(wildcards.sample, "biomart_mart"),
        dataset=lambda wildcards: metadata_value(wildcards.sample, "biomart_dataset"),
        id_filter=lambda wildcards: metadata_value(wildcards.sample, "biomart_id_filter"),
        id_attribute=lambda wildcards: metadata_value(wildcards.sample, "biomart_id_attribute"),
        strict="--strict" if config.get("biomart", {}).get("strict", False) else "",
    
    log:
        f"{LOGS}/{{sample}}/fetch_biomart.log"

    conda:
        "../envs/r.yaml"
    
    message:
        ""
        
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.annotations})" "$(dirname {log})"

        Rscript workflow/scripts/annotations/fetch_biomart_annotations.R \
          --input-gene-table {input.gene_table} \
          --sample {params.sample} \
          --host {params.host} \
          --mart {params.mart} \
          --dataset {params.dataset} \
          --id-filter {params.id_filter} \
          --id-attribute {params.id_attribute} \
          {params.strict} \
          --output {output.annotations} \
          > {log} 2>&1
        """
