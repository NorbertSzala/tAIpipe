"""
Snakemake workflow for constructing genome-specific cytosolic eukaryotic
ribosomal gene reference sets using KofamScan.

The workflow contains four sequential stages:

1. ``prepare_kofam_ribosomal_reference``
   Extracts cytosolic eukaryotic ribosomal KO identifiers from a local
   KEGG BRITE ko03011 hierarchy, validates the local KOfam resources and
   creates a restricted HAL profile list.

2. ``run_kofamscan_ribosome``
   Searches each genome-specific protein FASTA against the selected
   ribosomal KOfam HMM profiles.

3. ``parse_kofamscan_ribosome``
   Retains KO assignments passing predefined KOfam thresholds and maps
   protein identifiers to gene and CDS identifiers.

4. ``extract_ribosomal_reference_cds``
   Extracts nucleotide CDS records corresponding to significant ribosomal
   protein annotations.

The global ribosomal KO/HMM reference is created once and reused for all
genomes. KofamScan and downstream parsing are performed separately for each
sample.

The workflow assumes that protein identifiers are unique within each
protein FASTA and can be mapped unambiguously to CDS records using the
gene-protein mapping table.
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


# ---------------------------------------------------------------------------
# Output path patterns
# ---------------------------------------------------------------------------

KOFAM_DETAIL = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosome_detail.tsv"
)

KOFAM_HITS = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosome_significant_hits.tsv"
)

KOFAM_GENES = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosome_gene_annotations.tsv"
)

KOFAM_REFERENCE_GENE_IDS = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosomal_reference_gene_ids.txt"
)

KOFAM_REFERENCE_CDS_IDS = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosomal_reference_cds_ids.txt"
)

KOFAM_REFERENCE_CDS = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosomal_reference_cds.fna"
)

KOFAM_QC = (
    f"{PER_GENOME}/{{sample}}/kofamscan/"
    "ribosome_qc.tsv"
)

GENE_PROTEIN_MAP = (
    f"{PER_GENOME}/{{sample}}/tables/"
    "gene_protein_map.tsv"
)


def get_gene_protein_map(wildcards):
    """Return the expected gene-protein mapping path for a sample."""
    return GENE_PROTEIN_MAP.format(sample=wildcards.sample)


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
        detail=KOFAM_DETAIL,

    params:
        threshold_scale=KOFAM_THRESHOLD_SCALE,
        benchmark_dir=lambda wildcards: f"{BENCHMARKS}/{wildcards.sample}",

    threads:
        KOFAM_THREADS

    resources:
        mem_mb=KOFAM_MEM_MB

    log:
        f"{LOGS}/{{sample}}/kofamscan.log"

    benchmark:
        f"{BENCHMARKS}/{{sample}}/kofamscan.tsv"

    conda:
        "../envs/kofamscan.yaml"

    shell:
        """
        set -euo pipefail

        mkdir -p \
            "$(dirname {output.detail:q})" \
            "$(dirname {log:q})" \
            {params.benchmark_dir:q}

        tmpdir="$(
            mktemp -d \
            "${{TMPDIR:-/tmp}}/kofamscan.{wildcards.sample}.XXXXXX"
        )"

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
        detail=KOFAM_DETAIL,
        proteins=get_proteome,
        mapping=get_gene_protein_map,
        ko_table=RIBOSOME_KO_TABLE,

    output:
        hits=KOFAM_HITS,
        genes=KOFAM_GENES,
        gene_ids=KOFAM_REFERENCE_GENE_IDS,
        cds_ids=KOFAM_REFERENCE_CDS_IDS,
        qc=KOFAM_QC,

    params:
        gene_column=KOFAM["columns"].get(
            "gene_id",
            "gene_id",
        ),
        protein_column=KOFAM["columns"].get(
            "protein_id",
            "protein_id",
        ),
        cds_column=KOFAM["columns"].get(
            "cds_id",
            "cds_id",
        ),
        max_missing_proteome_count=config.get(
            "gene_protein_map",
            {},
        ).get("max_missing_proteome_count", 200),
        max_missing_proteome_fraction=config.get(
            "gene_protein_map",
            {},
        ).get("max_missing_proteome_fraction", 0.05),

    conda:
        "../envs/kofamscan.yaml"

    script:
        "../scripts/kofamscan/parse_kofamscan_ribosome.py"


rule extract_ribosomal_reference_cds:
    input:
        cds=get_cds,
        cds_ids=KOFAM_REFERENCE_CDS_IDS,

    output:
        cds=KOFAM_REFERENCE_CDS,

    conda:
        "../envs/kofamscan.yaml"

    script:
        "../scripts/kofamscan/extract_reference_cds.py"

rule select_top2_ribosomal_reference_cds:
    input:
        hits=f"{PER_GENOME}/{{sample}}/kofamscan/ribosome_significant_hits.tsv",
        cds=get_cds

    output:
        ids=f"{PER_GENOME}/{{sample}}/kofamscan/ribosomal_reference_cds_ids.top2_per_ko.txt",
        fasta=f"{PER_GENOME}/{{sample}}/kofamscan/ribosomal_reference_cds.top2_per_ko.fna"

    params:
        top_n=2,
        allow_empty=False

    log:
        stderr=f"{LOGS}/kofamscan/{{sample}}.select_top2_ribosomal_reference_cds.stderr.log"

    benchmark:
        f"{BENCHMARKS}/kofamscan/{{sample}}.select_top2_ribosomal_reference_cds.tsv"

    threads: 1

    resources:
        mem_mb=4000,
        disk_mb=20000,
        runtime=60

    conda:
        "../envs/python.yaml"

    message:
        "Selecting top 2 ribosomal CDS sequences per KO for {wildcards.sample}."

    script:
        "../scripts/kofamscan/select_top_n_reference_cds.py"