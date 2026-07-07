"""
Runs the R script that calculates per-gene codon-usage metrics for each genome. It combines validated CDS sequences with the corresponding tRNA copy-number profile and produces gene-level values such as ENC, CAI, FOP, tAI, GC, and GC3s, together with supporting codon-level files.
"""
rule codon_usage_metrics:
    input:
        cds=get_cds,
        aaa_count=f"{PER_GENOME}/{{sample}}/counted_codons/{{sample}}_aaa_counts.tsv",
        trna_qc=f"{PER_GENOME}/{{sample}}/qc/{{sample}}_trna_profile_qc.tsv",
        reference_cds = (
            f"{PER_GENOME}/{{sample}}/kofamscan/ribosomal_reference_cds.fna"
        )

    output:
        codon_counts=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_codon_counts.csv",
        enc=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_enc.csv",
        rscu=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_rscu.csv",
        cai=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_cai.csv",
        trna_weights=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_trna_weights.csv",
        tai=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_tai.csv",
        amino_acid_usage=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_amino_acid_usage.csv",
        fop=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_fop.csv",
        gc=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_gc.csv",
        gc3s=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_gc3s.csv",
        summary=f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv",
    
    params:
        domain=get_domain,
        sample=lambda wildcards: wildcards.sample,
        gcode=get_genetic_code,
        outdir=f"{PER_GENOME}/{{sample}}/codon_metrics",
    
    log:
        f"{LOGS}/{{sample}}/calculate_tai.log"

    conda:
        "../envs/r.yaml"

    message:
        "Counting codon usage metrics using R-Cubar package"

    shell:
        """
        set -euo pipefail 

        mkdir -p {params.outdir} $(dirname {log})
        
        Rscript workflow/scripts/metrics/calculate_tAI.R \
            --input {input.cds} \
            --reference_cds {input.reference_cds} \
            --outdir {params.outdir} \
            --genetic-code {params.gcode} \
            --trna {input.aaa_count} \
            --domain {params.domain} \
            --sample {params.sample} \
            > {log} 2>&1
        
        """