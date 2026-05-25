
rule codon_usage_metrics:
    """
    Count tAI.
    """

    input:
        cds = get_cds,
        aaa_count = f'{PER_GENOME}/{{sample}}/counted_codons/{{sample}}_aaa_counts.tsv', # trna counts - 


    
    output:
        summary = f'{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv'

    params:
        domain = get_domain,
        sample = lambda wildcards: wildcards.sample,
        gcode = get_genetic_code,
        outdir = f'{PER_GENOME}/{{sample}}/codon_metrics'
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
        
        Rscript workflow/scripts/calculate_tAI.R \
            --input {input.cds} \
            --outdir {params.outdir} \
            --genetic-code {params.gcode} \
            --trna {input.aaa_count} \
            --domain {params.domain} \
            --sample {params.sample} \
            > {log} 2>&1
        
        """