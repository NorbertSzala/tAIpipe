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

- Automated tRNA gene prediction using **tRNAscan-SE**
- CDS-based codon usage analysis, including **tAI**, **CAI**, **ENC**, **RSCU**, **GC**, and **GC3s** ([Brief definitions](https://www.genscript.com/gsfiles/tools/Index_Definition_of_GenRCA_Rare_Codon_Analysis_Tool.pdf))
- Support for genome-scale comparative analyses across organisms
- Designed for fungal datasets, but not restricted to fungi
- Snakemake-based workflow with reproducible per-genome outputs
- Outputs prepared for downstream plotting, statistical analysis, and dashboard integration

## ℹ️ Overview

`tAIpipe` is a bioinformatics workflow for studying codon usage and translational adaptation across genomic datasets. The pipeline combines genome/CDS input data, tRNA gene prediction, codon usage metrics, and structured output tables that can be used for comparative analyses and interactive visualization.

The main goal of this project is to integrate **tRNA Adaptation Index** analysis with organism-level metadata, such as taxonomy, genome size, lifestyle, and other biological features. This makes it possible to explore questions such as whether **tAI** depends on genome size, ecological niche, protein function, or the presence of selected domains.

The workflow is designed as a modular Snakemake pipeline. Each step produces intermediate files that can be inspected, reused, or extended in downstream analyses.

## What is **tAI**?

The **tRNA Adaptation Index (tAI)** measures how well codons in a coding sequence match the available pool of tRNAs in a cell. It uses tRNA gene copy numbers as a proxy for tRNA availability and accounts for codon–anticodon pairing, including wobble interactions.

Higher **tAI** values suggest that a gene uses codons recognized by more abundant tRNAs, which may support faster or more efficient translation. Lower **tAI** values may reflect locally slower translation, which can also be biologically meaningful, for example in protein folding or regulatory regions. **tAI** is therefore useful for studying codon optimization, translational efficiency, and evolutionary adaptation of coding sequences.

## ✍️ Authors

- [Norbert Szala](https://github.com/NorbertSzala)
- [Max Stróżyk](https://github.com/maxi7524)


## Quick start

```bash
git clone https://github.com/NorbertSzala/tAIpipe.git
cd tAIpipe
micromamba create -n snakemake -c conda-forge -c bioconda snakemake
micromamba activate snakemake
snakemake -n --profile workflow/profiles/test
snakemake --profile workflow/profiles/test
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
```


## 🚀 Usage

```bash
git clone https://github.com/NorbertSzala/tAIpipe.git
cd tAIpipe

snakemake --cores 4
````

To perform a dry run:

```bash
snakemake -n --cores 4 --printshellcmds
```

The main output files are generated per genome, for example:

```text
results/per_genome/<sample>/
├── trnascan/
├── counted_codons/
└── codon_usage/
```


### Execution backends

The workflow is designed to support two execution modes:

1. **Conda/Mamba mode** — current default and recommended mode on systems where Apptainer/Singularity is not available for unprivileged users.
2. **Apptainer/Singularity mode** — prepared as a reproducible containerized backend, but it requires a working system-level Apptainer/Singularity installation that can run containers without sudo.

Run with Conda:

```bash
snakemake --profile workflow/profiles/test
```

Run with Apptainer, if available:

```bash
snakemake --profile workflow/profiles/apptainer
```

## ⬇️ Installation

Clone the repository:

```bash
git clone https://github.com/NorbertSzala/tAIpipe.git
cd tAIpipe
```

Create the required environments. The pipeline uses both Python and R tools.

Basic requirements:

```text
Python
R
Snakemake
tRNAscan-SE
HMMER / hmmsearch
Biostrings
cubar
pandas
```

Example environment setup with micromamba:

```bash
micromamba create -n taipipe \
  -c conda-forge \
  -c bioconda \
  python \
  pandas \
  snakemake \
  trnascan-se \
  hmmer \
  r-base \
  r-readr \
  r-dplyr \
  r-optparse \
  bioconductor-biostrings \
  bioconductor-iranges \
  -y
```

Activate the environment:

```bash
micromamba activate taipipe
```

Install the R package `cubar` inside the environment:

```bash
Rscript -e 'install.packages("cubar", repos="https://cloud.r-project.org")'
```

## 🧩 Pipeline overview

The pipeline starts from a sample table defining organisms to analyse. Organisms can be specified by local genome/CDS files or by accession identifiers, depending on the configuration.

Main workflow steps:

1. Define organisms and input files in `config/config.yaml` and `data/metadata/samples.tsv`
2. Predict tRNA genes with `tRNAscan-SE`
3. Convert raw `tRNAscan-SE` output into a clean TSV table
4. Count amino acid–anticodon pairs, for example `Ala-CGC`
5. Calculate codon usage and adaptation metrics:

   * tAI
   * CAI
   * ENC
   * RSCU
   * GC
   * GC3s
6. Export tables for downstream plots, statistics, and dashboard integration

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

## 📊 Planned dashboard

The output tables are designed for integration with an interactive dashboard. Planned analyses include:

* tAI vs genome size
* tAI vs GC content
* tAI vs protein length
* tAI vs taxonomy
* tAI vs lifestyle or ecological niche
* tAI distributions across organisms
* codon usage differences between selected groups
* links between tAI and PFAM domain presence

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
