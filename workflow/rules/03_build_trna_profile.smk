"""
Transforms cleaned tRNAscan-SE predictions into the amino-acid and anticodon copy-number table required for tAI weight estimation. It also ensures that unsupported or non-elongator records are removed before the tRNA profile is passed to quality control and metric calculation.
"""
rule prepare_trna_codon_counts_to_tai:
    """
    Extract anticodon_id count data from cleaned tRNAscan-SE output

    Expected output:
        TSV with tRNA gene copy numbers in format:
        anticodon_id    count
        Ala-CGC         5

    """
    
    input:
        # tRNAscanse cleaned output
        clean_trnascan = f"{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.tsv"

    output:
        aaa_count = f'{PER_GENOME}/{{sample}}/counted_codons/{{sample}}_aaa_counts.tsv'


    log:
        f'{LOGS}/{{sample}}/prepare_trna_codon_counts_to_tai.log'

    conda:
        "../envs/python.yaml"
    
    message:
        "Preparing cleaned tRNAscan-SE data to interpretation in R script"

    shell:
        """
        mkdir -p $(dirname {output.aaa_count}) $(dirname {log})

        python3 workflow/scripts/preprocessing/prepare_trna_codon_counts_to_tai.py \
            -I {input.clean_trnascan} \
            -O {output.aaa_count} \
            > {log} 2>&1
        
        """