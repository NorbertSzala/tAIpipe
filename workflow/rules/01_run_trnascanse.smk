"""
Runs tRNAscan-SE independently for each genome. The rule applies the configured organism domain, genetic-code settings, sensitivity mode, and number of threads, and stores the raw prediction output for subsequent cleaning and counting.
"""

rule run_trnascanse:
    """
    Run tRNAscanSE-2.0 on given genome .fasta files. Predict tRNA genes in genome. Default option is to omit pseudo genes
    
    Intruction: https://gensoft.pasteur.fr/docs/trnascan/1.3.1/Manual.pdf

    Used flags:

    -E: eukaryotic mode (default)
    --max: maximum sensitivity mode - search using Infernat without hmm finding
    --output: path to main output
    --stats: path to output with stats
    --bed: path to output .bed file
    --gff: path to output .gff file
    --fasta: path to output .fasta file - saved predicted tRNA sequences
    --prefix: allows set {{sample}} as prefix
    --forceow: do not ask if overwrite
    --thread: set number of threads
    input.genome - input


    !IMPORTANT!

    in the output's column 'Anti codon' there are sequences of tRNA anticodons, f.e. AAT in DNA 5'-> 3'. They are not codons from coding sequence. 
    
    tRNAscan-SE reports tRNA-Ile with AAT anticodon, so biological-real-life sequence in tRNA is RNA 5'-AAU-3'. Proper codon in CDS is reverse complemented: DNA 5'-ATT-3' and that 
    """

    input:
        genome = get_genome
    
    output:
        out = f"{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.out",
        stats = f'{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.stats',
        bed = f'{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.bed',
        gff = f'{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.gff',
        fasta = f'{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.fasta',

    log:
        f'{LOGS}/{{sample}}/trnascanse.log'
    
    threads:
        config.get("trnascanse", {}).get("threads", 4)

    conda:
        "../envs/trnascanse.yaml"

    # container:
    #     CONTAINERS['trnascanse']

    params:
        sensitivity = lambda wildcards: "--max" if config.get("trnascanse", {}).get("max_sensitivity", True) else "",
        gencode_arg = get_trnascanse_code_arg,
        domain = get_trnascanse_domain_arg
        

    message:
        "Predicting tRNA-genes using tRNAscan-SE 2.0 tool."


    shell:
        """
        mkdir -p $(dirname {output.out}) $(dirname {log})  
        
        tRNAscan-SE \
            {params.domain} \
            {params.sensitivity} \
            {params.gencode_arg} \
            --output {output.out} \
            --stats {output.stats} \
            --bed {output.bed} \
            --gff {output.gff} \
            --fasta {output.fasta} \
            --prefix {wildcards.sample} \
            --forceow \
            --thread  {threads} \
            {input.genome} \
            > {log} 2>&1
        """
