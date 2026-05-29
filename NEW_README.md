# 🧬→🖥️→📊 tAIpipe

v0.4.0

Bioinformatics pipeline for tRNA adaptation index (**tAI**), codon usage, and tRNA gene analysis across genomic and CDS datasets, with outputs designed for downstream visualization and interactive dashboards.

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Snakemake](https://img.shields.io/badge/snakemake-039475?style=for-the-badge)
![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
![Power BI](https://img.shields.io/badge/power_bi-F2C811?style=for-the-badge&logo=powerbi&logoColor=black)
![Bioinformatics](https://img.shields.io/badge/bioinformatics-pipeline-green?style=for-the-badge)
![Status](https://img.shields.io/badge/status-in_development-orange?style=for-the-badge)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## 🌟 Highlights 

#TODO - insert here -> functionality / usage / For who
- ja bym tutaj opisał dokkładnie jaki jest cel tego (bo samo powiedzenie co jest nie jest potrzebne) 
- ogólnie jakie informacje możemy uzyskać i jaka jest ogólność tego

## ℹ️ Overview

#TODO 
- ja bym opisał tutaj dokładniej na czym polega pipeline oraz w jaki sposób można to zgeneralizować  
- czyli jakie dane (+- format) itd.
  - można już  tutaj rozbić na 
  - ### What is tAI
    - quick explanation of what we measure / detect 
  - ### How pipeline looks like  *this information from flow part
    - here we explain very briefly what is happening at each stage
  - #### More deeply into metric and scripts explanation 
    - breafly explanation per paragraph with link do deeper .md file 
  - #### what analysis includes (shortly and deeply analysis)
    - żeby tak dać - skrótowy opis i dać odnośni od tych dokladniejszych docsów
    - tutaj nas interesuje moment raport

### 📖 Further reading

Recommended papers and resources:

* [The tAI implementation concept](https://doi.org/10.1093/nar/gkg897)
* [Codon usage bias overview](https://doi.org/10.1007/s11033-021-06749-4)
* [Recent introduction to codon usage and translational regulation](https://doi.org/10.1038/s41467-024-52660-4)
* [Review on codon usage and translation dynamics](https://doi.org/10.1146/annurev-biophys-030722-020555)

## Usage
In this section we explain how to use our library.

#TODO - Norbert
- możemy zrobić mały tutoria w osobnym .md zeby właśnie na tym próbnym datasetcie to pokazać. 

### Environment configuration
In this sections we explain how to configure environment 

> UWAGA: Tutaj zamiast "Environment" można jakiś backend wstawić itd. 

#### Conda
`tAIpipe` leverages Snakemake's native environment virtualization. You do **not** need to manually install heavy software stacks (such as tRNAscan-SE, HMMER, or multi-package R environments) globally on your operating system. Instead, Snakemake dynamically reads individual rules and isolates software stacks into local, dedicated project directories inside the repository root path under: `resources/conda/<environment-hash>/`.


To configure your base orchestrator environment properly using a high-performance Mamba framework, execute the following commands:

```bash
# 1. Initialize a fully compliant Conda/Mamba core deployment (Miniforge)
# Skip this if a full Conda core environment manager is already configured
curl -L -O "[https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh](https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh)"
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
source ~/miniforge3/bin/activate
conda init bash

# Restart your terminal session after running the initialization script
## Case: linux
soruce ~/.bashrc 
## Case: general
## Close and open terminal 

# 2. Compile a minimal environment container dedicated strictly to Snakemake
mamba create -n snakemake_env -c conda-forge snakemake -y
mamba activate snakemake_env
```

> ⚠️ **CRITICAL ARCHITECTURAL NOTE:** Snakemake 9.x strictly requires a fully standard-compliant Conda/Mamba metadata interface (such as Miniforge) to evaluate execution prefixes and dependency maps. Using a standalone, stripped-down lightweight deployment (like pure global `micromamba`) will cause metadata parsing crashes during runtime, specifically throwing `KeyError: 'conda_prefix'` or exit status 127.


#TODOTERAZ 

#### Apptainer
#TODO - wyłumacz jak to ustawić bo z tego nie korzystałem i nie wiem jak to działa. Potem przeczytam i sprawdzę czy działa 

#### Apptainer - test
Apptainer (formerly Singularity) allows you to execute the pipeline inside fully containerized environments instead of dynamically compiling Conda packages. This is highly useful for HPC clusters or multi-user systems where package downloading or internet access is restricted during runtime.

To utilize Apptainer, ensure that the system-level Apptainer binary is installed and that your user can run unprivileged containers (without `sudo`). 

You can execute the entire workflow using the containerized backend profile:
```bash
snakemake --profile workflow/profiles/apptainer
```

When this backend is invoked, Snakemake completely skips local Conda setups and pulls pre-built Docker/Singularity images defined inside the rule configurations, mounting your workspace automatically inside the container instance.

### Snakemake configuration
#TODO - Norbert:
- tutaj wytłumaczyć w jaki sposób ustawić configi w snakemake, i zwrócić uwage na przybliżone zużycie zasobów obliczeniowych, żeby użytkownik +- wiedział ile RAMu i czasu potrzebuje na wykonanie tego. Nawet jak to jest mało to trzeba skomentować
- jak mamy te dane testowe, to właśnie żeby pokazać na nich, ALBO najlepiej zrobić osobny .md z szybkim tutorialem jak to na tej mniejszej próbce uruchomic

### Workflow execution
After configuring your environment profiles and verifying your input configuration files, you can proceed to trigger the processing pipeline.

#### Standard Pipeline Execution

This is the default command to start actual data processing. Snakemake reads the test profile configuration, evaluates the execution DAG (Directed Acyclic Graph), maps missing target outputs, spawns up to 8 parallel worker threads, and launches the tasks sequentially:

```bash
snakemake --profile workflow/profiles/test
```

#### Passive Dry Run Verification

Before committing heavy compute resources, it is highly recommended to inspect the execution plan without producing any files or running shell jobs. Adding the -n (or --dry-run) flag forces Snakemake to display exactly which rules will be triggered, which input files are bound, and the explicit reasoning behind each job execution (e.g., missing outputs or updated upstream files):

```bash
snakemake -n --profile workflow/profiles/test
```

## 📁 Scripts and rules

### Scripts

```text
scripts/
├── calculate_tAI.R
│   └── Calculates tAI, CAI, ENC, RSCU, GC, GC3s, and related codon usage metrics.
│
├── convert_trnascanse_output_to_tsv.py
│   └── Converts raw tRNAscan-SE output into a clean TSV table.
│
├── downsample_genomes.sh
│   └── Extracts regions around selected tRNA genes to speed up testing and prototyping.
│
├── extract_first_contig.sh
│   └── Extracts the first contig from genome FASTA files for lightweight test runs.
│
└── prepare_trna_codon_counts_to_tai.py
    └── Converts cleaned tRNAscan-SE output into amino acid–anticodon copy-number counts.
        Example output format:
            anticodon_id    count
            Ala-CGC         4
            Gly-TCC         7
            Leu-CAG         3
```

### Snakemake rules

```text
workflow/rules/
├── trnascan_rule.smk
│   └── Runs tRNAscan-SE and predicts tRNA genes.
│
├── clean_tRNAscanSE_output.smk
│   └── Cleans and reformats tRNAscan-SE output.
│
├── prepare_trna_codon_counts_to_tai_rule.smk
│   └── Counts amino acid–anticodon pairs for tAI calculation.
│
└── codon_usage_metrics_rule.smk
    └── Calculates tAI and additional codon usage metrics.
```







## 💭 Feedback and contributing

Suggestions, bug reports, and feature requests are welcome.

Open an issue if you find:

* incorrect parsing of input files
* unsupported genome/CDS formats
* problems with tRNAscan-SE output conversion
* missing codon usage metrics
* ideas for new dashboard features

Contributions are welcome through pull requests or private message.



## ✍️ Authors

- [Norbert Szala](https://github.com/NorbertSzala)
- [Max Stróżyk](https://github.com/maxi7524)


#TODO - tutaj poźniej rozpisć kto za co był odpowiedzialny 