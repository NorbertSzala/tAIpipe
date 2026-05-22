# Workflow overview

`tAIpipe` is a Snakemake workflow for comparative tRNA adaptation and codon usage analysis.

The workflow is designed to process multiple fungal genomes and their corresponding CDS FASTA files. It predicts tRNA genes, prepares amino acid–anticodon copy-number tables, and calculates tAI and codon usage metrics.

---

## Main steps

1. Read workflow configuration from `config/config.yaml`
2. Validate the main configuration using `workflow/schemas/config.schema.yaml`
3. Read sample metadata from `config/samples.tsv`
4. Filter samples using the `include` column
5. Resolve genome and CDS input files using filename patterns from the sample table
6. Run `tRNAscan-SE` on genome or downsampled genome FASTA
7. Convert raw `tRNAscan-SE` output into clean TSV
8. Remove or keep pseudogenes depending on config settings
9. Remove undetermined tRNAs if requested
10. Count amino acid–anticodon pairs
11. Calculate tAI and codon usage metrics
12. Export per-genome and aggregated result tables

---

## Current workflow structure

```text
workflow/
├── Snakefile
├── rules/
│   ├── trnascan_rule.smk
│   ├── clean_tRNAscanSE_output.smk
│   ├── prepare_trna_codon_counts_to_tai_rule.smk
│   └── codon_usage_metrics_rule.smk
├── envs/
│   ├── python.yaml
│   ├── r-cubar.yaml
│   └── trnascanse.yaml
├── schemas/
│   └── config.schema.yaml
├── profiles/
│   ├── default/
│   └── test/
└── containers/
    ├── containers.yaml
    └── defs/
```

---

## Core rules

| rule                               | purpose                                                                  |
| ---------------------------------- | ------------------------------------------------------------------------ |
| `run_trnascanse`                   | predicts tRNA genes from genome FASTA using `tRNAscan-SE`                |
| `clean_tRNAscanSE_output`          | converts raw `tRNAscan-SE` output into a clean TSV table                 |
| `prepare_trna_codon_counts_to_tai` | counts amino acid–anticodon pairs from cleaned tRNA predictions          |
| `codon_usage_metrics`              | calculates tAI, CAI, ENC, RSCU, GC, GC3s and related codon usage metrics |

---

## Main input files

The workflow uses two main input sources:

### 1. Configuration file

```text
config/config.yaml
```

This file defines paths and analysis parameters, for example:

```yaml
trnascanse:
  threads: 8
  max_sensitivity: true
  domain: "-E"

trnascanse_clean:
  keep_pseudo: false
  remove_nnn: true
```

### 2. Sample table

The sample table is defined in `config/config.yaml` under:

```yaml
paths:
  samples: ...
```

Expected columns include:

```text
sample
species
genetic_code
accession
domain
kingdom
phylum
lifestyle
genome
cds
include
```

The `include` column controls whether a sample is used in the workflow. Only rows with:

```text
include == True
```

are selected.

---

## Path resolution

Genome and CDS files are not hardcoded directly in rules. Instead, the workflow reads filename patterns from the sample table and resolves the matching files with helper functions.

Example:

```python
def get_genome(wildcards):
    pattern = samples_df.loc[wildcards.sample, "genome"]
    return resolve_single_file(f"{DATA_GENOME}/{pattern}")
```

The same approach is used for CDS files:

```python
def get_cds(wildcards):
    pattern = samples_df.loc[wildcards.sample, "cds"]
    return resolve_single_file(f"{DATA_CDS}/{pattern}")
```

This allows the sample table to contain patterns such as:

```text
GCA_000219625.1*.fna
```

The workflow requires that each pattern resolves to exactly one file. If no file or multiple files are found, the workflow stops with an error.

---

## Genetic code handling

The genetic code is sample-specific and is read from the `genetic_code` column in the sample table.

Example helper:

```python
def get_genetic_code(wildcards):
    return int(samples_df.loc[wildcards.sample, "genetic_code"])
```

This allows different organisms to use different genetic codes, for example:

| organism                      | genetic code |
| ----------------------------- | -----------: |
| most fungi                    |            1 |
| `Schizosaccharomyces pombe`   |            3 |
| `Candida albicans`-like cases |           12 |

---

## tRNAscan-SE execution

The `run_trnascanse` rule runs `tRNAscan-SE` on the genome FASTA.

Important parameters are read from `config/config.yaml`:

```yaml
trnascanse:
  threads: 8
  max_sensitivity: true
  domain: "-E"
```

The rule uses:

| flag        | meaning                  |
| ----------- | ------------------------ |
| `-E`        | eukaryotic mode          |
| `--max`     | maximum sensitivity mode |
| `--output`  | main tRNAscan-SE output  |
| `--stats`   | statistics output        |
| `--bed`     | BED output               |
| `--gff`     | GFF output               |
| `--fasta`   | predicted tRNA sequences |
| `--prefix`  | sample-specific prefix   |
| `-gencode`  | genetic code             |
| `--forceow` | overwrite existing files |
| `--thread`  | number of threads        |

Important biological note:

The `Anticodon` column in the `tRNAscan-SE` output contains tRNA anticodons, not CDS codons. For example, a reported DNA anticodon `AAT` corresponds to the RNA anticodon `AAU`; the matching CDS codon is the reverse complement, `ATT`.

---

## Amino acid–anticodon counts

After cleaning the `tRNAscan-SE` output, the workflow counts amino acid–anticodon pairs.

The output has the format:

```text
anticodon_id    count
Ala-CGC         5
Arg-ACG         2
Gly-GCC         4
```

This format is intended for downstream tAI calculation with R/cubar.

Pseudogenes and undetermined tRNAs can be removed according to:

```yaml
trnascanse_clean:
  keep_pseudo: false
  remove_nnn: true
```

---

## Test dataset generation

A lightweight test dataset can be created from selected NCBI accession numbers using:

```bash
./workflow/scripts/create_test_dataset.sh \
  -C data/CDS/ \
  -G data/genome/ \
  -A "GCA_000219625.1,GCA_000006255.1,GCA_000002945.2,GCA_000149205.2,GCA_000181695.1,GCA_002104945.1,GCA_016906535.1,GCA_905067625.1"
```

The script searches for files matching:

```text
data/CDS/{accession}*.fna
data/genome/{accession}*.fna
```

For each accession, it creates:

```text
resources/test_data/
├── CDS/
│   └── ACCESSION_first_contig_CDS.fna
└── genome/
    └── ACCESSION_first_contig_1Mbp.fna
```

Genome downsampling strategy:

1. Take the first genome contig.
2. Truncate it to the first 1,000,000 bp.
3. Try to select CDS records matching the first contig ID.
4. If no matching CDS records are found, use a fallback subset of the first CDS records.

This produces a small but realistic dataset for fast workflow testing.

---

## Metadata for test dataset

A reduced metadata table can be generated from the full `Main_Dataset.tsv` by keeping only rows corresponding to selected accessions.

Recommended pattern:

```bash
mkdir -p resources/test_data/metadata

{
  head -n 1 resources/metadata/Main_Dataset.tsv

  for accession in \
    GCA_000219625.1 \
    GCA_000006255.1 \
    GCA_000002945.2 \
    GCA_000149205.2 \
    GCA_000181695.1 \
    GCA_002104945.1 \
    GCA_016906535.1 \
    GCA_905067625.1
  do
    grep "$accession" resources/metadata/Main_Dataset.tsv
  done
} > resources/test_data/metadata/test_metadata.tsv
```

This preserves the header and writes a clean test metadata file.

---

## Execution environments

The current working execution mode is:

```text
Snakemake + Conda/Mamba environments
```

The workflow keeps tool-specific Conda environment files in:

```text
workflow/envs/
```

Typical environments:

| environment       | purpose                                                 |
| ----------------- | ------------------------------------------------------- |
| `python.yaml`     | Python scripts, e.g. pandas-based preprocessing         |
| `trnascanse.yaml` | `tRNAscan-SE` execution                                 |
| `r-cubar.yaml`    | R, Biostrings, IRanges, cubar and plotting dependencies |

The recommended execution model is:

1. Activate only the controller environment containing Snakemake.
2. Let Snakemake create and activate rule-specific environments.

Example:

```bash
micromamba activate project
snakemake --profile workflow/profiles/test
```

Do not manually switch between `trnascanse`, `r-cubar`, and Python environments while running the workflow. Those should be managed by Snakemake through the `conda:` directive.

---

## Conda/Mamba mode

Current practical execution mode:

```bash
snakemake --profile workflow/profiles/test
```

or explicitly:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --profile workflow/profiles/test \
  --use-conda \
  --printshellcmds
```

The profile should include:

```yaml
use-conda: true
conda-frontend: mamba
conda-prefix: resources/conda
```

The directory:

```text
resources/conda/
```

is used for Snakemake-created environments and should not be committed to Git.

Recommended `.gitignore` entries:

```gitignore
.snakemake/
resources/conda/
resources/containers/
results/
```

---

## Apptainer/Singularity container support

The workflow has initial support for Apptainer/Singularity containers through:

```text
workflow/containers/
```

Suggested structure:

```text
workflow/containers/
├── containers.yaml
└── defs/
    ├── python.def
    └── r-cubar.def
```

Purpose:

| path                                  | purpose                                      |
| ------------------------------------- | -------------------------------------------- |
| `workflow/containers/containers.yaml` | stores image names or local `.sif` paths     |
| `workflow/containers/defs/`           | stores Apptainer definition files            |
| `resources/containers/`               | stores generated or downloaded `.sif` images |

The container definitions are kept for reproducibility and future deployment. However, on the current server, unprivileged Apptainer execution is not available, producing errors related to user namespaces and permissions. Therefore, the active execution backend is Conda/Mamba.

Important note:

Apptainer/Singularity containers require a working system-level installation that can run containers without `sudo`. Snakemake does not run jobs through `sudo apptainer`. If this command fails:

```bash
apptainer exec image.sif echo OK
```

then Snakemake with `--use-apptainer` will also fail.

Current recommendation:

```text
Use Conda/Mamba backend on this server.
Keep Apptainer/Singularity definitions as a documented reproducibility backend.
```

---

## Minimal dry run

For the test profile:

```bash
snakemake -n --profile workflow/profiles/test
```

Expected DAG for the current test dataset:

```text
8 × run_trnascanse
8 × convert_trnascanse_output_to_tsv
8 × prepare_trna_codon_counts_to_tai
1 × all
= 25 jobs
```

---

## Minimal real run

```bash
snakemake --profile workflow/profiles/test
```

For easier debugging, first run with one core:

```bash
snakemake --profile workflow/profiles/test --cores 1 --printshellcmds
```

Then run with more cores:

```bash
snakemake --profile workflow/profiles/test --cores 4 --printshellcmds
```
