"""
Defines the rules that construct the canonical analysis tables used by downstream statistics and plotting steps. It coordinates the creation of gene-level, genome-level, and codon-level tables from per-genome metrics, metadata, annotations, and quality-control outputs.
"""

import shlex

ANNOTATION_TABLE = config.get("paths", {}).get("external_annotations", "")
GO_ENRICHMENT_ENABLED = config.get("go_enrichment", {}).get("enabled", False)
REQUIRE_GO_TERMS = config.get("go_enrichment", {}).get("require_go_terms", GO_ENRICHMENT_ENABLED)


def optional_external_annotations(wildcards):
    return [ANNOTATION_TABLE] if ANNOTATION_TABLE else []


def external_annotation_argument(wildcards):
    if not ANNOTATION_TABLE:
        return ""

    return (
        "--annotation-table "
        f"{shlex.quote(ANNOTATION_TABLE)}"
    )


def require_go_argument(wildcards):
    return "--require-go-terms" if REQUIRE_GO_TERMS else ""

"""Build the canonical gene-level feature table."""
rule build_gene_features:
    input:
        summaries=expand(
            f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv",
            sample=SAMPLES,
        ),
        trna_qc=expand(
            f"{PER_GENOME}/{{sample}}/qc/{{sample}}_trna_profile_qc.tsv",
            sample=SAMPLES,
        ),
        metadata_dataset=config["paths"]["metadata_dataset"],
        external_annotations=optional_external_annotations,

    output:
        table=config["paths"]["gene_features"],

    params:
        per_genome_dir=PER_GENOME,
        cds_dir=DATA_CDS,
        annotation_arg=external_annotation_argument,
        require_go_arg=require_go_argument,

    log:
        f"{LOGS}/tables/build_gene_features.log"

    conda:
        "../envs/r.yaml"

    shell:
        """
        set -euo pipefail

        mkdir -p \
            "$(dirname {output.table:q})" \
            "$(dirname {log:q})"

        Rscript workflow/scripts/tables/build_gene_features.R \
            --metadata-dataset {input.metadata_dataset:q} \
            --per-genome-dir {params.per_genome_dir:q} \
            --cds-dir {params.cds_dir:q} \
            {params.annotation_arg} \
            {params.require_go_arg} \
            --output {output.table:q} \
            > {log:q} 2>&1
        """


# ----------------------------------------------------------------

rule build_genome_summary:
    input:
        gene_features=config["paths"]["gene_features"],
        trna_qc=expand(
            f"{PER_GENOME}/{{sample}}/qc/{{sample}}_trna_profile_qc.tsv",
            sample=SAMPLES,
        ),
        genomes=[
            resolve_single_file(f"{DATA_GENOME}/{samples_df.loc[sample, 'genome']}")
            for sample in SAMPLES
        ],
        metadata_dataset=config["paths"]["metadata_dataset"],

    output:
        table=config["paths"]["genome_summary"],

    params:
        genome_dir=DATA_GENOME,
        per_genome_dir=PER_GENOME,

    log:
        f"{LOGS}/tables/build_genome_summary.log"

    conda:
        "../envs/r.yaml"

    message:
        "Building the canonical genome-level summary table"

    shell:
        """
        set -euo pipefail

        mkdir -p \
            "$(dirname {output.table:q})" \
            "$(dirname {log:q})"

        Rscript workflow/scripts/tables/build_genome_summary.R \
            --gene-features {input.gene_features:q} \
            --metadata-dataset {input.metadata_dataset:q} \
            --genome-dir {params.genome_dir:q} \
            --per-genome-dir {params.per_genome_dir:q} \
            --output {output.table:q} \
            > {log:q} 2>&1
        """

# ----------------------------------------------------------------

rule build_codon_profiles:
    input:
        codon_counts=expand(
            f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_codon_counts.csv",
            sample=SAMPLES,
        ),
        rscu=expand(
            f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_rscu.csv",
            sample=SAMPLES,
        ),
        reference_rscu=expand(
            f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_reference_rscu.csv",
            sample=SAMPLES,
        ),
        trna_weights=expand(
            f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_trna_weights.csv",
            sample=SAMPLES,
        ),
        metadata_dataset=config["paths"]["metadata_dataset"],

    output:
        table=config["paths"]["codon_profiles"],

    params:
        per_genome_dir=PER_GENOME,

    log:
        f"{LOGS}/tables/build_codon_profiles.log"

    conda:
        "../envs/r.yaml"

    message:
        "Building the canonical genome-by-codon profile table"

    shell:
        """
        set -euo pipefail

        mkdir -p \
            "$(dirname {output.table:q})" \
            "$(dirname {log:q})"

        Rscript workflow/scripts/tables/build_codon_profiles.R \
            --metadata-dataset {input.metadata_dataset:q} \
            --per-genome-dir {params.per_genome_dir:q} \
            --output {output.table:q} \
            > {log:q} 2>&1
        """
