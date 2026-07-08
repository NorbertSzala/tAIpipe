"""
Defines quality-control checks performed before and after codon-metric calculation. It validates the completeness and structure of each tRNA profile and can also verify that metric outputs contain the expected columns, identifiers, and sufficient finite values.
"""
rule validate_trna_profile:
    input:
        counts=f"{PER_GENOME}/{{sample}}/counted_codons/{{sample}}_aaa_counts.tsv",

    output:
        qc=f"{PER_GENOME}/{{sample}}/qc/{{sample}}_trna_profile_qc.tsv",

    params:
        sample=lambda wildcards: wildcards.sample,
        mode=config["trna_qc"]["mode"],
        min_total_trnas=config["trna_qc"]["min_total_trnas"],
        min_unique_anticodons=config["trna_qc"]["min_unique_anticodons"],
        min_amino_acids=config["trna_qc"]["min_amino_acids"],

    log:
        f"{LOGS}/{{sample}}/validate_trna_profile.log"

    conda:
        "../envs/python.yaml"

    message:
        "Validating tRNA profile for {wildcards.sample}"

    shell:
        """
        set -euo pipefail

        mkdir -p \
            "$(dirname {output.qc:q})" \
            "$(dirname {log:q})"

        python workflow/scripts/qc/validate_trna_profile.py \
            --input {input.counts:q} \
            --output {output.qc:q} \
            --sample {params.sample:q} \
            --mode {params.mode:q} \
            --min-total-trnas {params.min_total_trnas} \
            --min-unique-anticodons {params.min_unique_anticodons} \
            --min-amino-acids {params.min_amino_acids} \
            > {log:q} 2>&1
        """
