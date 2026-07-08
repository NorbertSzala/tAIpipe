"""
Snakemake workflow for constructing genome-specific cytosolic eukaryotic
ribosomal gene reference sets using KofamScan.

The workflow contains four sequential stages:

1. ``prepare_kofam_ribosomal_reference``
   Extracts cytosolic eukaryotic ribosomal KO identifiers from a local
   KEGG BRITE ko03011 hierarchy, validates the local KOfam resources and
   creates a restricted HAL profile list.

2. ``run_kofamscan_ribosome``
   Searches each genome-specific protein FASTA against only the selected
   ribosomal KOfam HMM profiles.

3. ``parse_kofamscan_ribosome``
   Retains KO assignments passing predefined KOfam thresholds, maps protein
   identifiers to gene identifiers and generates gene-level annotations and
   quality-control metrics.

4. ``extract_ribosomal_reference_cds``
   Extracts nucleotide coding sequences for the selected ribosomal genes,
   producing the reference FASTA used in downstream codon-usage and CAI
   calculations.

The global ribosomal KO/HMM reference is created once and reused for every
genome. KofamScan and downstream parsing are performed separately for each
sample.

The workflow assumes that protein identifiers are unique within every
protein FASTA and can be mapped unambiguously to gene identifiers.
"""

KOFAM = config["kofamscan"]

KOFAM_KO_LIST = KOFAM["database"]["ko_list"]
KOFAM_PROFILES = KOFAM["database"]["profiles_dir"]

RIBOSOME_BRITE = KOFAM["ribosome_reference"]["brite_json"]
RIBOSOME_KO_TABLE = KOFAM["ribosome_reference"]["ko_table"]
RIBOSOME_HAL = KOFAM["ribosome_reference"]["hal"]

KOFAM_THREADS = int(
    KOFAM["execution"].get("threads", 8)
)
KOFAM_MEM_MB = int(
    KOFAM["execution"].get("mem_mb", 16000)
)
KOFAM_THRESHOLD_SCALE = float(
    KOFAM["execution"].get("threshold_scale", 1.0)
)
KOFAM_DETAIL = (
    f"{PER_GENOME}/{{sample}}/kofamscan/ribosome_detail.tsv"
)

KOFAM_HITS = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosome_significant_hits.tsv"
)

KOFAM_GENES = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosome_gene_annotations.tsv"
)

KOFAM_REFERENCE_IDS = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosomal_reference_cds_ids.txt"
)

KOFAM_REFERENCE_CDS = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosomal_reference_cds.fna"
)

def format_sample_path(pattern, wildcards):
    return pattern.format(sample=wildcards.sample)


def get_proteins(wildcards):
    return get_proteome(wildcards)


def get_cds(wildcards):
    return format_sample_path(
        KOFAM["inputs"]["cds"],
        wildcards,
    )


def get_gene_protein_map(wildcards):
    return (
        f"{PER_GENOME}/{wildcards.sample}/tables/gene_protein_map.tsv"
    )


rule prepare_kofam_ribosomal_reference:
    input:
        brite_json=RIBOSOME_BRITE,
        ko_list=KOFAM_KO_LIST,

    output:
        ko_table=RIBOSOME_KO_TABLE,
        hal=RIBOSOME_HAL,

    params:
        profiles_dir=KOFAM_PROFILES,
        require_thresholds=KOFAM[
            "ribosome_reference"
        ].get("require_thresholds", True),

    conda:
        "../envs/kofamscan.yaml"

    script:
        "../scripts/kofamscan/prepare_ribosomal_reference.py"


rule run_kofamscan_ribosome:
    input:
        proteins=get_proteome,
        profile=RIBOSOME_HAL,
        ko_list=KOFAM_KO_LIST,

    output:
        detail=(
            "results/per_genome/{sample}/kofamscan/"
            "ribosome_detail.tsv"
        ),
    
    params:
        threshold_scale = KOFAM_THRESHOLD_SCALE,
        
    threads:
        KOFAM_THREADS

    resources:
        mem_mb=KOFAM_MEM_MB

    log:
        f"{LOGS}/{{sample}}/kofamscan.log"

    benchmark:
        f"{BENCHMARKS}/{{sample}}//kofamscan.tsv"

    conda:
        "../envs/kofamscan.yaml"
    
    shell:
        """
        set -euo pipefail

        mkdir -p \
            "$(dirname {output.detail:q})" \
            "$(dirname {log:q})" \
            "$(dirname {benchmark:q})"

        tmpdir="$(mktemp -d "${{TMPDIR:-/tmp}}/kofamscan.XXXXXX")"

        cleanup() {{
            rm -rf "$tmpdir"
        }}

        trap cleanup EXIT

        exec_annotation \
            --profile {input.profile:q} \
            --ko-list {input.ko_list:q} \
            --cpu {threads} \
            --tmp-dir "$tmpdir" \
            --format detail-tsv \
            --threshold-scale {params.threshold_scale} \
            --no-report-unannotated \
            -o {output.detail:q} \
            {input.proteins:q} \
            > {log:q} 2>&1

        test -s {output.detail:q}
        """


rule parse_kofamscan_ribosome:
    input:
        detail=(
            "results/per_genome/{sample}/kofamscan/"
            "ribosome_detail.tsv"
        ),
        proteins=get_proteins,
        mapping=get_gene_protein_map,
        ko_table=RIBOSOME_KO_TABLE,

    output:
        hits=(
            "results/per_genome/{sample}/kofamscan/"
            "ribosome_significant_hits.tsv"
        ),

        genes=(
            "results/per_genome/{sample}/kofamscan/"
            "ribosome_gene_annotations.tsv"
        ),

        gene_ids=(
            "results/per_genome/{sample}/kofamscan/"
            "ribosomal_reference_gene_ids.txt"
        ),

        qc=(
            "results/per_genome/{sample}/kofamscan/"
            "ribosome_qc.tsv"
        ),

    params:
        gene_column=KOFAM["columns"].get(
            "gene_id",
            "gene_id",
        ),

        protein_column=KOFAM["columns"].get(
            "protein_id",
            "protein_id",
        ),

    conda:
        "../envs/kofamscan.yaml"

    script:
        "../scripts/kofamscan/parse_kofamscan_ribosome.py"


rule extract_ribosomal_reference_cds:
    input:
        cds=get_cds,
        cds_ids=(
            f"{PER_GENOME}/{{sample}}/kofamscan/ribosomal_reference_cds_ids.txt")

    output:
        cds=(
            "{PER_GENOME}/{{sample}}/kofamscan/"
            "ribosomal_reference_cds.fna"
        ),

    conda:
        "../envs/kofamscan.yaml"
        
    script:
        "../scripts/kofamscan/extract_reference_cds.py"
