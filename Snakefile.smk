# Load config file from yaml file
configfile: "config.yaml"

# Helper variables for easier acces
DATA_CDS = config["paths"]["data_cds"]
DATA_DOWN = config["paths"]["data_downsampled"]
DATA_GENOME = config["paths"]["data_genome"]
DATA_MAINTABLE = config["paths"]["data_main_table"]
DATA_PLOTS = config["paths"]["data_plots"]
DATA_TABLES = config["paths"]["data_tables"]
REPORT_FILE = ["paths"]["report_file"]


# --- Main Rule ---
rule all:
    input:
        config["paths"]["report_file"]

# ----------------------------------
# ---  Process data to count tAI ---
# ----------------------------------
rule run_trnascanse:
    """
    Run tRNAscanSE-2.0 on given genome .fasta files. Predict tRNA genes in genome.
    """

    input:
    
    output:

    shell:


rule codonM:
    """
    Count codon statistics
    """
    input:
        # tRNAscanse output

    output:
    
    shell:


rule codonW:
    """
    Count codons frequencies
    """

    input:
        #cds
    
    output:

    shell:


rule count_and_transform_trnas:
    '''
    Program reads table from tRNAscan-SE output and count previously complemented codons. Output is in format needed in gtAI (R).

    !!!!!!!!!!! CHECK THIS INFORMATION BELOW!!!!!!!!!!!!!!!
    Note: input Anti Codon is in DNA (as seq in input tRNAscan-SE)
    '''
    input:
        # trnascanse ouptut
    
    output:

    shell:


rule count_tAI:
    """
    Count tAI.
    """

    input:
        # trnascanse ouptut, codonw output, codonm output, count_codons.py output
    
    output:

    shell: