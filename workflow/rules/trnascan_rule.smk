rule run_trnascanse:
    """
    Run tRNAscanSE-2.0 on given genome .fasta files. Predict tRNA genes in genome. Default option is to omit pseudo genes
    
    Intruction: https://gensoft.pasteur.fr/docs/trnascan/1.3.1/Manual.pdf

    Used flags:

    -E: eukariotic mode (default)
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


    #TODO IMPORTANT: Later, add flag: --max


    !IMPORTANT!

    in the output's column 'Anti codon' there are sequences of tRNA anticodons, f.e. AAT in DNA 5'-> 3'. They are not codons from coding sequence. 
    
    tRNAscan-SE reports tRNA-Ile with AAT anticodon, so biological-real-life sequence in tRNA is RNA 5'-AAU-3'. Proper codon in CDS is reverse complemented: DNA 5'-ATT-3' and that 
    """

    input:
        genome = get_first_contig
    
    output:
        out = f"{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.out",
        stats = f'{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.stats',
        bed = f'{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.bed',
        gff = f'{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.gff',
        fasta = f'{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.fasta',

    log:
        f'{LOGS}/{{sample}}/trnascan.log'

    threads:
        config.get("trnascanse_threads", 4)

    conda:
        "../../envs/trnascan.yaml"

    shell:
        """
        mkdir -p $(dirname {output.out}) $(dirname {log})  
        
        tRNAscan-SE \
            -E \
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


# TODO: dodaj obsluge alternate genetic code - flaga --gencode