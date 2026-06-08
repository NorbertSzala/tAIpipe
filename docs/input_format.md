# Input format

This document describes the core input files required to run `tAIpipe`. It defines the expected FASTA inputs, sample table structure, metadata fields, and file-matching rules used by the Snakemake workflow.

For a general workflow overview, see [`docs/workflow_overview.md`](workflow_overview.md).
For metrics description, see [`docs/metrics.md`](metrics.md)

---

## 1. Input architecture

`tAIpipe` requires three main input components:

1. **Genome FASTA files**
   Used by `tRNAscan-SE` to predict tRNA genes.

2. **CDS FASTA files**
   Used to calculate codon usage metrics, including tAI, CAI, ENC, RSCU, GC, and GC3s.

3. **Sample table**
   Defines sample identifiers, input file paths or filename patterns, genetic codes, taxonomy, lifestyle categories, and inclusion flags.

The workflow uses separate configuration files for test and production runs:

| Run mode   | Config file               | Sample table              | Input data                  |
| ---------- | ------------------------- | ------------------------- | --------------------------- |
| Test       | `config/config_test.yaml` | `config/samples_test.tsv` | `resources/test_data/`      |
| Production | `config/config.yaml`      | `config/samples.tsv`      | `data/genome/`, `data/CDS/` |

---

## 2. Genome FASTA files

Genome FASTA files are required because tRNA genes are non-coding RNA genes and are not represented in CDS FASTA files. They are used as input for `tRNAscan-SE`.

Typical locations:

```text
resources/test_data/genome/
data/genome/
```

Example file:

```text
head GCA_000002945.2_first_contig_1Mbp.fna
 >contig_or_chromosome_id
ATGCGT...
```

---

## 3. CDS FASTA files

CDS FASTA files are used to calculate codon usage and translational adaptation metrics.

Typical locations:

```text
resources/test_data/CDS/
data/CDS/
```

Example file:

```text
head GCA_000002945.2_first_contig_CDS.fna
>gene_or_cds_id
ATGGCT...
```

The sequences should be valid coding sequences compatible with the genetic code specified in the sample table.

---

## 4. Sample table

The sample table is the main control file used by the workflow. It defines which samples are processed and where their genome/CDS files are located.

Test sample table:

```text
config/samples_test.tsv
```

Only rows with:

```text
include = True
```

are processed.

---

## 5. Required sample table columns

| Column         | Description                                                                                                             |
| -------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Short unique sample identifier used in output paths. Example: `Spombe`.                                                 |
| `species`      | Species name. Example: `Schizosaccharomyces_pombe`.                                                                     |
| `genetic_code` | NCBI genetic code ID used for codon usage and tAI-related calculations. Example: `1`, `12`.                             |
| `domain`       | Biological domain used to select the appropriate `tRNAscan-SE` mode. Expected values: `eukarya`, `bacteria`, `archaea`. |
| `accession`    | Genome assembly accession ID. Example: `GCA_000002945.2`.                                                               |
| `kingdom`      | Broad taxonomic group. Example: `Fungi`.                                                                                |
| `phylum`       | Phylum-level taxonomy. Example: `Ascomycota`.                                                                           |
| `lifestyle`    | Ecological or biological lifestyle category used for downstream grouping.                                               |
| `genome`       | Genome FASTA filename, path, or pattern.                                                                                |
| `cds`          | CDS FASTA filename, path, or pattern.                                                                                   |
| `include`      | Boolean flag controlling whether the sample should be processed. Expected values: `True` or `False`.                    |

---

## 6. Optional sample table columns

The following columns may be present but are not required by the core workflow:

| Column              | Description                                                                                |
| ------------------- | ------------------------------------------------------------------------------------------ |
| `bed`               | Optional BED file with genomic coordinates. Reserved for downstream or auxiliary analyses. |
| `downsampled_fasta` | Optional reference to a reduced FASTA file used during development or testing.             |
| `first_contig`      | Optional reference to the first-contig FASTA file used for lightweight test datasets.      |

These columns are useful for documenting how test data were generated, but they are not required for the main tAI/codon usage workflow.

---

## 7. Example sample table

Example `config/samples_test.tsv`:

```tsv
sample	species	genetic_code	domain	accession	kingdom	phylum	lifestyle	genome	cds	include
Spombe	Schizosaccharomyces_pombe	1	eukarya	GCA_000002945.2	Fungi	Ascomycota	nectar_tap_saprotroph	resources/test_data/genome/GCA_000002945.2_first_contig_1Mbp.fna	resources/test_data/CDS/GCA_000002945.2_first_contig_CDS.fna	True
Calbicans	Candida_albicans	12	eukarya	GCA_000182965.3	Fungi	Ascomycota	pathogen	resources/test_data/genome/GCA_000182965.3_first_contig_1Mbp.fna	resources/test_data/CDS/GCA_000182965.3_first_contig_CDS.fna	True
```

---

---

## 8. Matching constraints

The workflow intentionally uses strict file matching.

- If a pattern does not match any file, Snakemake raises a `MissingInputException`.
- 
  ```text
  GCA_999999999.1*.fna
  ```


- If a pattern matches more than one file, the workflow should stop with an explicit error.

  ```text
  GCA_0002*.fna
  ```

Each `genome` and `cds` entry should resolve to one unique FASTA file.

---

## 9. Genetic codes

The `genetic_code` column stores the NCBI genetic code ID used for each sample.

For standard nuclear code:

```text
1
```

For selected non-standard codes, `tAIpipe` uses tRNAscan-SE-compatible genetic-code files from:

```text
data/genetic_codes/trnascanse/
```

Example files:

```text
12.gcode
3.gcode
4.gcode
6.gcode
```

For NCBI genetic code `1`, no additional tRNAscan-SE genetic-code file is required.

The .gcode file does not need to contain the complete genetic code. Instead, it should specify only the codons that differ from the standard genetic code. See the tRNAscan-SE documentation for details.

---

## 10. Metadata used during aggregation and reports

The sample table contains both execution-level fields and biological grouping fields.

Execution-level fields include:

```text
sample
genetic_code
domain
genome
cds
include
```

Biological grouping fields include:

```text
species
accession
kingdom
phylum
lifestyle
```

The core computational rules use execution-level fields to build the Snakemake DAG. Biological grouping fields are mainly used during aggregation, plotting, and report generation.

Optional external metadata tables may be added later if additional biological variables are needed, such as:

```text
class
order
family
genome_size
microenvironment
host_association
```

Such metadata should be merged during the aggregation/reporting stage using stable identifiers such as `sample`, `accession`, or `species`.

---

## 11. Test data

The repository includes a small test dataset:

```text
resources/test_data/
├── CDS/
├── genome/
└── metadata/
```

This dataset is used by:

```text
config/config_test.yaml
config/samples_test.tsv
```

It is intended for quick workflow validation:

```bash
snakemake -n --profile workflow/profiles/test
snakemake --profile workflow/profiles/test
```

---

## 12. Production data

Production runs use:

```text
config/config.yaml
config/samples.tsv
```

and usually read full genome/CDS files from:

```text
data/genome/
data/CDS/
```

Before running production mode, check that all paths in `config/samples.tsv` resolve to existing files.

Recommended dry run:

```bash
snakemake -n --profile workflow/profiles/production
```

---

## 13. Common input problems

| Problem                        | Likely cause                      | Fix                                       |
| ------------------------------ | --------------------------------- | ----------------------------------------- |
| `MissingInputException`        | Genome or CDS path does not exist | Check `genome` and `cds` columns          |
| Multiple files matched         | Pattern is too broad              | Use a more specific filename or full path |
| tRNAscan-SE genetic-code error | Missing `.gcode` file             | Check `data/genetic_codes/trnascanse/`    |
| Sample not processed           | `include` is `False`              | Set `include` to `True`                   |
| Empty output table             | CDS file is empty or invalid      | Check CDS FASTA content                   |
| Translation/codon warnings     | Wrong genetic code                | Verify `genetic_code` for the sample      |
