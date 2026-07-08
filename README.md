# 🧬→🖥️→📊 tAIpipe

v0.6.0

`tAIpipe` is a Snakemake workflow for comparative analysis of tRNA gene repertoires, codon usage and translational adaptation. For every included genome, it predicts tRNA genes, constructs a cytosolic ribosomal reference set, calculates gene-level codon metrics and produces canonical tables for downstream statistics.

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Snakemake](https://img.shields.io/badge/snakemake-039475?style=for-the-badge)
![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
![Power BI](https://img.shields.io/badge/power_bi-F2C811?style=for-the-badge&logo=powerbi&logoColor=black)
![Bioinformatics](https://img.shields.io/badge/bioinformatics-pipeline-green?style=for-the-badge)
![Status](https://img.shields.io/badge/status-in_development-orange?style=for-the-badge)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## 🌟 Highlights

- Automated tRNA gene prediction using **tRNAscan-SE**
- CDS-based codon usage analysis, including **tAI**, **CAI**, **ENC**, **RSCU**, **GC**, and **GC3s** ([Brief definitions](https://www.genscript.com/gsfiles/tools/Index_Definition_of_GenRCA_Rare_Codon_Analysis_Tool.pdf))
- Support for genome-scale comparative analyses across organisms
- Designed for fungal datasets, but not restricted to fungi
- Snakemake-based workflow with reproducible per-genome outputs
- Outputs prepared for downstream plotting, statistical analysis, and dashboard integration

## ℹ️ Overview

`tAIpipe` is a modular Snakemake workflow for studying codon usage and translational adaptation across genomic datasets. The pipeline combines genome/CDS input data, tRNA gene prediction, codon usage metrics, and structured output tables that can be used for comparative analyses and visualization.

The main goal of this project is to integrate **tRNA Adaptation Index** analysis with organism-level metadata, such as taxonomy, genome size, lifestyle, and other biological features. This makes it possible to explore whether **tAI** depends on genome size, ecological niche, protein function, or selected genomic features.

Detailed documentation is split into separate files:

| Guide                                                    | Contents                                     |
| -------------------------------------------------------- | -------------------------------------------- |
| [`docs/getting_started.md`](docs/getting_started.md)     | Installation, configuration and first run    |
| [`docs/input_format.md`](docs/input_format.md)           | Required inputs and FASTA/header constraints |
| [`docs/workflow_overview.md`](docs/workflow_overview.md) | Rule order and data flow                     |
| [`docs/script_reference.md`](docs/script_reference.md)   | Input and output of every active script      |
| [`docs/output_format.md`](docs/output_format.md)         | Output directory tree and table schemas      |
| [`docs/metrics.md`](docs/metrics.md)                     | Metric and statistical interpretation        |

## What is **tAI**?

The **tRNA Adaptation Index (tAI)** measures how well codons in a coding sequence match the available pool of tRNAs in a cell. It uses tRNA gene copy numbers as a proxy for tRNA availability and accounts for codon–anticodon pairing, including wobble interactions.

Higher **tAI** values suggest that a gene uses codons recognized by more abundant tRNAs, which may support faster or more efficient translation. Lower **tAI** values may reflect locally slower translation, which can also be biologically meaningful, for example in protein folding or regulatory regions.

## ✍️ Authors

- [Norbert Szala](https://github.com/NorbertSzala)
- [Max Stróżyk](https://github.com/maxi7524)



## Requirements

- Snakemake
- Conda, Mamba or Micromamba
- Genome, CDS and protein FASTA files for every sample
- Local KOfam `ko_list` and profile HMM database
- KEGG BRITE `ko03011` JSON used to define cytosolic eukaryotic ribosomal proteins

Rule-specific software is installed from `workflow/envs/`. The R environment installs the `cubar` package through `workflow/envs/r.post-deploy.sh`.

A working Conda/Mamba/Micromamba installation is recommended.

## Quick start

### 0. Download

```bash
git clone https://github.com/NorbertSzala/tAIpipe.git
cd tAIpipe
```

### 1. Validate the test workflow

```bash
snakemake -n --profile workflow/profiles/test
```

### 2. Run the test workflow

```bash
snakemake --profile workflow/profiles/test
```

### 3. Run the production workflow

Edit `config/config.yaml` and `config/samples.tsv` (adjust to your paths, samples), then run:

```bash
snakemake -n --profile workflow/profiles/production
snakemake --profile workflow/profiles/production
```

A direct invocation without a profile is also supported:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --use-conda \
  --cores 16 \
  --rerun-incomplete \
  --printshellcmds
```

## Information flow

```text
Genome FASTA
   ↓
tRNAscan-SE
   ↓
clean TSV
   ↓
amino acid–anticodon counts
   ↓
R/cubar metrics
   ↓
per-genome summary tables
   ↓
aggregated metrics table
```

## 🚀 Usage

### Test dry run

```bash
snakemake -n --profile workflow/profiles/test
```

### Test run

```bash
snakemake --profile workflow/profiles/test
```


### Production run

```bash
snakemake --profile workflow/profiles/production
```

### Manual run without profile

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --cores 8 \
  --use-conda \
  --printshellcmds \
  --rerun-incomplete
```


## Execution backends

The workflow supports two execution modes:

1. **Conda/Mamba mode** — default and recommended for local/HPC use.
2. **Apptainer/Singularity mode** — optional containerized backend for systems with Apptainer/Singularity support. NOTE: THIS FEATURE IS STILL UNDER DEVELOPMENT

Run with Conda:

```bash
snakemake --profile workflow/profiles/test
```

NOTE: THIS FEATURE IS STILL UNDER DEVELOPMENT
Run with Apptainer, if available:

```bash
snakemake --profile workflow/profiles/apptainer
```

## 🧩 Pipeline overview

Main workflow steps:

1. Define organisms and input files in `config/config.yaml` and `config/samples.tsv`.
2. Predict tRNA genes with `tRNAscan-SE`.
3. Convert raw `tRNAscan-SE` output into a clean TSV table.
4. Count amino acid–anticodon pairs, for example `Ala-CGC`.
5. Calculate codon usage and adaptation metrics.
6. Aggregate per-genome summaries into final tables.
7. Export tables for downstream plots, statistics, and dashboard integration.
8. 

```mermaid
flowchart TD
    %% =========================================
    %% Per-sample genomic analysis workflow
    %% =========================================

    subgraph SAMPLE["Per-sample workflow"]
        direction TB

        A["Genome FASTA"]
            --> B["1. tRNAscan-SE"]

        B --> C["Raw tRNAscan-SE output"]

        C --> D["2. Clean tRNAscan-SE output"]

        D --> E["Clean tRNA table"]

        E --> F["3. Count amino acid–anticodon pairs"]

        F --> G["tRNA anticodon count table"]

        H["CDS FASTA"]

        G --> I["4. Codon usage metrics<br/>and tAI calculation"]
        H --> I

        I --> J["Per-sample codon usage summary"]

        J --> K["5. Aggregate metrics"]
    end

    K --> L["Final multi-sample tables and plots"]

    %% =========================================
    %% Styling
    %% =========================================

    classDef input fill:#0f172a,stroke:#334155,color:#ffffff,stroke-width:1px;
    classDef process fill:#1e3a8a,stroke:#1d4ed8,color:#ffffff,stroke-width:1px;
    classDef intermediate fill:#e2e8f0,stroke:#94a3b8,color:#0f172a,stroke-width:1px;
    classDef output fill:#065f46,stroke:#10b981,color:#ffffff,stroke-width:1px;

    class A,H input;
    class B,D,F,I,K process;
    class C,E,G,J intermediate;
    class L output;
```


## 📁 Scripts and rules

Main Snakemake rules are stored in:

```text
workflow/rules/
```

Main analysis scripts are stored in:

```text
workflow/scripts/
```

Important components:

```text
workflow/rules/01_run_trnascanse.smk
workflow/rules/02_clean_trnascanse.smk
workflow/rules/03_build_trna_profile.smk
workflow/rules/07_compute_codon_metrics.smk
workflow/rules/aggregate_metrics.smk
workflow/rules/plots.smk
```

```text
workflow/scripts/calculate_tAI.R
workflow/scripts/convert_trnascanse_output_to_tsv.py
workflow/scripts/prepare_trna_codon_counts_to_tai.py
workflow/scripts/aggregate_metrics.py
```

## 📊 Outputs

The main output files are generated per genome and then aggregated.

```text
results/
├── per_genome/
│   └── <sample>/
│       ├── trnascan/
│       ├── trna_counts/
│       └── codon_metrics/
└── aggregated/
    ├── tables/
    ├── plots/
    └── reports/
```

Main final table:

```text
results/aggregated/tables/all_genomes_metrics.tsv
```

### Main outputs

| File                                        | Unit                        | Purpose                                          |
| ------------------------------------------- | --------------------------- | ------------------------------------------------ |
| `results/tables/gene_features.tsv`          | one CDS/gene                | Gene metadata, QC and codon-usage metrics        |
| `results/tables/genome_summary.tsv`         | one genome                  | Genome-level summaries and QC                    |
| `results/tables/codon_profiles.tsv`         | one genome × codon          | Codon counts, RSCU, CAI weights and tRNA weights |
| `results/statistics/gene_feature_tests.tsv` | one binary feature          | Mixed-model association with within-genome tAI   |
| `results/statistics/genome_group_tests.tsv` | one metric × group variable | Genome-level Wilcoxon or Kruskal–Wallis tests    |
| `results/statistics/go_enrichment.tsv`      | one GO term × tAI tail      | Genome-stratified GO enrichment                  |


See `docs/output_format.md` for detailed output descriptions.


## 📖 Further reading

Recommended papers and resources:

* [The tAI implementation concept](https://doi.org/10.1093/nar/gkg897)
* [Codon usage bias overview](https://doi.org/10.1007/s11033-021-06749-4)
* [Recent introduction to codon usage and translational regulation](https://doi.org/10.1038/s41467-024-52660-4)
* [Review on codon usage and translation dynamics](https://doi.org/10.1146/annurev-biophys-030722-020555)

## 💭 Feedback and contributing

Suggestions, bug reports, and feature requests are welcome.

Open an issue if you find:

* incorrect parsing of input files
* unsupported genome/CDS formats
* problems with tRNAscan-SE output conversion
* missing codon usage metrics
* ideas for new dashboard features

Contributions are welcome through pull requests or private message.


## Participation

#### Norbert Szala

- Designed the main concept and technical structure of the workflow.
- Implemented the Snakemake workflow organization, including:
  - configuration files in `config/*.yaml`,
  - execution profiles in `workflow/profiles/`,
  - input validation schemas in `workflow/schemas/`,
  - Conda environment definitions and lock files in `workflow/envs/locks/`,
  - repository directory structure and output organization.
- Implemented reproducibility-oriented workflow components, including Conda-based execution and initial Apptainer/Singularity container support.
- Developed core data-processing scripts for converting raw inputs into codon usage and tAI-related metrics, including:
  - parsing and cleaning `tRNAscan-SE` outputs,
  - preparing amino acid–anticodon count tables,
  - calculating codon usage metrics with R and the `cubar` package,
- Integrated the main analysis steps into the Snakemake pipeline.
- Prepared and maintained the main README and project documentation in `docs/`.
- Prepared test-oriented workflow configuration and small example-data execution mode.
- All plots inside `workflow/scripts` structure
- Prepare scripts counting all used metrics, including preparation to calculate CAI (and choosing ribosomal reference set using KEGG database)

#### Max Stróżyk
- Implementation of `.html` layer with plots inside

## License

GPL-3.0. See `LICENSE`.
