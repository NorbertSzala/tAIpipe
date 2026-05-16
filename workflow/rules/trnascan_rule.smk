rule run_trnascanse:
    """
    Run tRNAscanSE-2.0 on given genome .fasta files. Predict tRNA genes in genome. Default option is to omit pseudo genes
    
    Intruction: https://gensoft.pasteur.fr/docs/trnascan/1.3.1/Manual.pdf
    """

    input:
        genome = get_genome
    
    output:
        out = f"{PER_GENOME}/{{sample}}/trnascan/trnascan.out",
        stats = f'{PER_GENOME}/{{sample}}/trnascan/trnascan.stats'
    
    log:
        f'{LOGS}/{{sample}}/trnascan.log'

    conda:
        "../../envs/trnascan.yaml"

    shell:
        """
        mkdir -p $(dirname {output.out}) $(dirname {log})  
        
        tRNAscan-SE \
            -E \
            -o {output.out} \
            -m {output.stats} \
            {input.genome} \
            > {log} 2>&1
        """


