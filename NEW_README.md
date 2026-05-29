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


***

## ℹ️ Overview

#TODO 
- ja bym opisał tutaj dokładniej na czym polega pipeline oraz w jaki sposób można to zgeneralizować  
- czyli jakie dane (+- format) itd.
  - można już  tutaj rozbić na 



### What is tAI
The **tRNA Adaptation Index (tAI)** is a biosubstitutive metric used to estimate the translational efficiency of protein-coding genes. 

Instead of relying on costly global quantification of cellular tRNA molecules, tAI utilizes tRNA gene copy numbers extracted from genomic sequences as a approximation for tRNA availability. By integrating these copy numbers with domain-specific codon-anticodon wobble pairing. 

Higher **tAI** values suggest that a gene uses codons recognized by more abundant tRNAs, which may support faster or more efficient translation. Lower **tAI** values may reflect locally slower translation, which can also be biologically meaningful, for example in protein folding or regulatory regions. **tAI** is therefore useful for studying codon optimization, translational efficiency, and evolutionary adaptation of coding sequences.

### What the Analysis Includes
<!-- `tAIpipe` offers both a high-level overview and a deep-dive comprehensive analysis of codon preferences:
- **Global Profiling:** Quantifying total tRNA gene pools across various fungal genomes.
- **Translational Adaptation (tAI):** Scoring individual gene sequences based on translational optimization.
- **Codon Usage Bias (CUB):** Calculating classic evolutionary metrics such as the Frequency of Optimal Codons (FOP), Effective Number of Codons (ENC), and Relative Synonymous Codon Usage (RSCU).
- **Compositional Bias:** Assessing GC and synonymous third-position GC content (GC3s) to untangle mutational bias from natural selection.
- **Macro-Evolutionary Insights:** Aggregating per-genome metrics with ecological and taxonomic metadata (e.g., lifestyle, phylum) to uncover adaptive patterns across the fungal kingdom.

#TODO 
- tutaj napisać to w takim "uogólnionym" formacie, to jest że dla dowolnego gatunku możemy to ustawić 
- oraz bardziej napisać w formacie co robimy (skrótowo ale konkretnie) 
- 
- 
- 
- -->


### Pipeline explanation

#### Introduction to Pipeline Architecture
`tAIpipe` follows a strict **per-sample** design pattern. Processing samples individually ensures parallel computational efficiency and localized error handling before multi-sample data consolidation.

To navigate the workflow outputs and configurations, we define the following abstract path structures governed by `config/config.yaml`:
- `{DATA_GENOME}`: Directory containing raw or downsampled whole-genome FASTA files (`.fna`).
- `{DATA_CDS}`: Directory containing protein-coding sequence FASTA files (`.fna`).
- `{PER_GENOME}`: Target root directory for isolated, sample-specific results.
- `{AGGREGATED}`: Final directory where individual metrics are compiled, plotted, and summarized.

**Core Inputs Breakdown:**
- **Genome Files (`{DATA_GENOME}`)**: Required because tRNA genes are non-coding RNAs scattered throughout intergenic and intronic regions. We need the full chromosomal context to find them.
- **CDS Files (`{DATA_CDS}`)**: Required to extract the exact frequency of the 61 sense codons within the translated part of the genome.
- **Sample Table (`data/tutorial_data/input/metadata/samples_test.tsv`)**: The file containing experiment configuration for Snakemake, defining sample taxonomy, required genetic codes, file-naming patterns, and execution flags (`include`).
- **Metadata Master Table (`data/tutorial_data/input/metadata/test_dataset.tsv`)**: A comprehensive phenotypic matrix containing downstream categorical factors (lifestyles, microenvironments) used exclusively during final data aggregation.


#### Pipeline graph
```text
{DATA_GENOME} (FASTA)              {DATA_CDS} (FASTA)
                 │                                │
                 ▼                                │
      [ 1. tRNAscan-SE ]                          │
                 │                                │
                 ▼                                │
         Raw Text Output                          │
                 │                                │
                 ▼                                │
    [ 2. clean_tRNAscanSE_output ]                │
                 │                                │
                 ▼                                │
          Cleaned TSV Table                       │
                 │                                │
                 ▼                                │
 [ 3. prepare_trna_codon_counts_to_tai ]          │
                 │                                │
                 ▼                                │
      Anticodon Count Table                       │
                 │                                │
                 └───────────────┬────────────────┘
                                 │
                                 ▼
                    [ 4. codon_usage_metrics ]
                                 │
                                 ▼
                     {PER_GENOME}/<sample>/
                     (tAI, FOP, GC, Summary TSV)
                                 │
                                 ▼
                   [ 5. aggregate_and_report ]
                                 │
                                 ▼
                       {AGGREGATED}/reports/
                       (Interactive Dashboard & .md)
```

For a detailed map of the repository directories, code components, and environment configurations, see the [Project Map Explanation](docs/workflow_overview.md).


### 📖 Further reading

Recommended papers and resources:

* [The tAI implementation concept](https://doi.org/10.1093/nar/gkg897)
* [Codon usage bias overview](https://doi.org/10.1007/s11033-021-06749-4)
* [Recent introduction to codon usage and translational regulation](https://doi.org/10.1038/s41467-024-52660-4)
* [Review on codon usage and translation dynamics](https://doi.org/10.1146/annurev-biophys-030722-020555)

## Usage
In this section we explain how to use our library.

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

#### Apptainer
#TODO - wyłumacz jak to ustawić bo z tego nie korzystałem i nie wiem jak to działa. Potem przeczytam i sprawdzę czy działa 

<!-- #### Apptainer - test

#TODO - u mnie poszła taka propozycja, ale nie wiem jak to sprawdzić i co to dokładnie jest 
Apptainer (formerly Singularity) allows you to execute the pipeline inside fully containerized environments instead of dynamically compiling Conda packages. This is highly useful for HPC clusters or multi-user systems where package downloading or internet access is restricted during runtime.

To utilize Apptainer, ensure that the system-level Apptainer binary is installed and that your user can run unprivileged containers (without `sudo`). 

You can execute the entire workflow using the containerized backend profile:
```bash
snakemake --profile workflow/profiles/apptainer
```

When this backend is invoked, Snakemake completely skips local Conda setups and pulls pre-built Docker/Singularity images defined inside the rule configurations, mounting your workspace automatically inside the container instance. -->

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

<!-- 
DEPRACATED - przeniosłem całość do docs/workflow_overview.md, ale możemy zrobić osobny skrypt który będzie to przechowywał i wtedy to tutaj wyświetlić. 
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
``` -->


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