rule prepare_trna_codon_counts_to_tai:
    """
    Extract anticodon_id count data from cleaned tRNAscan-SE output
    """
    
    input:
        # tRNAscanse cleaned output
        clean_trnascan = f"{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.tsv"

    output:
        aaa_count = f'{PER_GENOME}/{{sample}}/counted_codons/{{sample}}_aaa_counts.tsv'


    log:
        f'{LOGS}/{{sample}}/prepare_trna_codon_counts_to_tai.log'

    conda:
        "workflow/envs/python.yaml"
    
    shell:
        """
        mkdir -p $(dirname {ou  tput.aaa_count}) $(dirname {log})

        python3 workflow/scripts/prepare_trna_codon_counts_to_tai.py \
            -I {input.clean_trnascan} \
            -O {output.aaa_count} \
            > {log} 2>&1
        
        """