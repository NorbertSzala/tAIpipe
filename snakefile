# Load config file from yaml file
configfile: "./config/config.yaml"

import pandas as pd
from glob import glob

# ----------------------
# --- Path variables ---
# ----------------------
# Helper path variables for easier acces

DATA_CDS = config["paths"]["data_cds"]
DATA_DOWN = config["paths"]["data_downsampled"]
DATA_GENOME = config["paths"]["data_genome"]

PER_GENOME = config['paths']['per_genome']
LOGS = config['paths']['logs']

DATA_MAINTABLE = config["paths"]["main_dataset"]
DATA_PLOTS = config["paths"]["plots"]
DATA_TABLES = config["paths"]["tables"]
REPORT_FILE = config["paths"]["reports"]


# ---------------------
# --- Samples table ---
# ---------------------
# Read samples from tsv
samples_df = pd.read_csv(config['paths']['samples'], sep = '\t')
samples_df = samples_df.set_index("sample")
SAMPLES = samples_df.index.tolist()


# -------------------------------
# ------- Helper functions ------
# -------------------------------
def resolve_single_file(pattern):
    matches = glob(pattern) # get list of path to proper file based on pattern

    if len(matches) == 0:
        raise FileNotFoundError(f"No file found for pattern: {pattern}")
    
    if len(matches)>1:
        raise ValueError(f'Multiple files found for pattern: {pattern}')
    
    return matches[0]


def get_genome(wildcards):
    pattern = samples_df.loc[wildcards.sample,'genome']
    return resolve_single_file(f"{DATA_GENOME}/{pattern}")
    

def get_cds(wildcards):
    pattern = samples_df.loc[wildcards.sample, 'cds']
    return resolve_single_file(f'{DATA_CDS}/{pattern}')


def get_downsampled_fasta(wildcards):
    pattern = samples_df.loc[wildcards.sample]['downsampled_fasta']    
    return resolve_single_file(f'{DATA_DOWN}/{pattern}')

def get_first_contig(wildcards):
    pattern = samples_df.loc[wildcards.sample]['first_contig']    
    return resolve_single_file(f'{DATA_DOWN}/{pattern}')


def get_bed(wildcards):
    pattern = samples_df.loc[wildcards.sample]['bed']
    return resolve_single_file(f'{DATA_DOWN}/{pattern}')




# ----------------------
# ---  Include rules ---
# ----------------------
include: "workflow/rules/trnascan_rule.smk"


# -----------------
# --- Main Rule ---
# -----------------

rule all:
    input:
        expand(
            f'{PER_GENOME}/{{sample}}/trnascan/trnascan.out', sample=SAMPLES
        )



# ----------------------------------
# ---  Process data to count tAI ---
# ----------------------------------






# rule count_and_transform_trnas:
#     '''
#     Program reads table from tRNAscan-SE output and count previously complemented codons. Output is in format needed in gtAI (R).

#     !!!!!!!!!!! CHECK THIS INFORMATION BELOW!!!!!!!!!!!!!!!
#     Note: input Anti Codon is in DNA (as seq in input tRNAscan-SE)
#     '''
#     input:
#         # trnascanse ouptut
    
#     output:

#     shell:


