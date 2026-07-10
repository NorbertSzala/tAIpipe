"""
Runs the R script that calculates per-gene codon-usage metrics for each genome. It combines validated CDS sequences with the corresponding tRNA copy-number profile and produces gene-level values such as ENC, CAI, FOP, tAI, GC, and GC3s, together with supporting codon-level files.
"""
rule codon_usage_metrics:
    input:
        cds=get_cds,
        aaa_count=f"{PER_GENOME}/{{sample}}/counted_codons/{{sample}}_aaa_counts.tsv",
        trna_qc=f"{PER_GENOME}/{{sample}}/qc/{{sample}}_trna_profile_qc.tsv",
        reference_cds=lambda wildcards: (
            f"{PER_GENOME}/{wildcards.sample}/kofamscan/"
            f"ribosomal_reference_cds.{REFERENCE_SUFFIX}.fna"
        )

    output:
        codon_counts=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_codon_counts.csv",
        enc=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_enc.csv",
        rscu=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_rscu.csv",
        reference_rscu = f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_reference_rscu.csv",
        cai=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_cai.csv",
        trna_weights=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_trna_weights.csv",
        tai=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_tai.csv",
        amino_acid_usage=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_amino_acid_usage.csv",
        fop=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_fop.csv",
        gc=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_gc.csv",
        gc3s=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_gc3s.csv",
        summary=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv",
        metric_qc=f"{PER_GENOME}/{{sample}}/qc/{{sample}}_metric_qc.tsv",
    
    params:
        domain=get_domain,
        sample=lambda wildcards: wildcards.sample,
        gcode=get_genetic_code,
        outdir=f"{PER_GENOME}/{{sample}}/codon_metrics",
        min_reference_cds = config.get('ribosomal_reference', {}).get("min_reference_genes", 20),
        min_finite_tai_fraction=config.get("trna_qc", {}).get("min_finite_tai_fraction", 0.90),
        min_used_codon_coverage=config.get("trna_qc", {}).get("min_used_codon_coverage", 0.95),
        qc_mode=config.get("trna_qc", {}).get("mode", "warn"),
    
    log:
        f"{LOGS}/{{sample}}/calculate_tai.log"

    conda:
        "../envs/r.yaml"

    message:
        "Counting codon usage metrics using R-Cubar package"

    shell:
        """
        set -euo pipefail 

        mkdir -p {params.outdir} $(dirname {output.metric_qc}) $(dirname {log})
        
        Rscript workflow/scripts/metrics/calculate_tAI.R \
            --input {input.cds} \
            --reference-cds {input.reference_cds} \
            --outdir {params.outdir} \
            --genetic-code {params.gcode} \
            --trna {input.aaa_count} \
            --domain {params.domain} \
            --sample {params.sample} \
            --min-reference-cds {params.min_reference_cds} \
            --qc-output {output.metric_qc} \
            --min-finite-tai-fraction {params.min_finite_tai_fraction} \
            --min-used-codon-coverage {params.min_used_codon_coverage} \
            --qc-mode {params.qc_mode} \
            > {log} 2>&1
        
        """