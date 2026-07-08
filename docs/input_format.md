# Input format

This document describes the core input files required to run `tAIpipe`. It defines the expected FASTA inputs, sample table structure, metadata fields, and file-matching rules used by the Snakemake workflow.

For a general workflow overview, see [`docs/workflow_overview.md`](workflow_overview.md).
For metrics description, see [`docs/metrics.md`](metrics.md)

---

## 1. Input architecture

Each included sample requires three sequence files:

| Input         | Used for                           | Required properties                                                                   |
| ------------- | ---------------------------------- | ------------------------------------------------------------------------------------- |
| Genome FASTA  | tRNA prediction and genome summary | Non-empty nucleotide FASTA; unique sequence IDs recommended                           |
| CDS FASTA     | Codon metrics and gene metadata    | Nucleotide CDS records; unique IDs; NCBI-style header attributes strongly recommended |
| Protein FASTA | KofamScan                          | Protein record IDs must be unique and match CDS `[protein_id=...]` values             |

The workflow uses separate configuration files for test and production runs:

| Run mode   | Config file               | Sample table              | Input data                  |
| ---------- | ------------------------- | ------------------------- | --------------------------- |
| Test       | `config/config_test.yaml` | `config/samples_test.tsv` | `resources/test_data/`      |
| Production | `config/config.yaml`      | `config/samples.tsv`      | `data/genome/`, `data/CDS/` |

### Sample table

The sample table path is configured as:

```yaml
paths:
  metadata_dataset: config/samples.tsv
```

Required columns:

| Column         | Meaning                                                                                       |
| -------------- | --------------------------------------------------------------------------------------------- |
| `sample`       | Unique output-safe identifier matching `^[A-Za-z0-9_.-]+$`                                    |
| `species`      | Species label copied to final tables                                                          |
| `genetic_code` | NCBI translation table ID used by `cubar`; also selects an optional tRNAscan-SE `.gcode` file |
| `accession`    | Assembly or dataset accession                                                                 |
| `domain`       | `eukarya`, `eukaryota`, `bacteria` or `archaea`; used for tRNA wobble weights                 |
| `kingdom`      | Taxonomic metadata                                                                            |
| `phylum`       | Default genome-level grouping variable                                                        |
| `lifestyle`    | Default genome-level grouping variable                                                        |
| `genome`       | Filename or glob pattern relative to `paths.data_genome`                                      |
| `cds`          | Filename or glob pattern relative to `paths.data_cds`                                         |
| `proteome`     | Filename or glob pattern relative to `paths.data_proteome`                                    |
| `include`      | Boolean-like value controlling sample inclusion                                               |


#### File matching

Every pattern must match exactly one file:

```text
0 matches  -> workflow stops
1 match    -> file is used
>1 matches -> workflow stops
```

Prefer accession-specific patterns such as:

```text
GCA_000002945.2*.fna
```
---

## 2. FASTA identifier requirements

### CDS FASTA

The mapping script expects NCBI-style attributes in CDS descriptions:

```text
>lcl|... [gene=abc1] [locus_tag=SPBC1] [protein_id=NP_000001.1] [protein=...]
```

Identifier roles:

| Identifier          | Source                                                  | Use                                      |
| ------------------- | ------------------------------------------------------- | ---------------------------------------- |
| `seq_id` / `cds_id` | First FASTA header token (`SeqRecord.id`)               | Gene-level metric key and CDS extraction |
| `protein_id`        | `[protein_id=...]`                                      | Links CDS to protein FASTA and KofamScan |
| `gene_id`           | `[locus_tag=...]`, then `[gene=...]`, then `protein_id` | Canonical biological gene label          |

CDS records without `[protein_id=...]` are skipped by `build_gene_protein_map.py`. Because strict matching is enabled, every mapped protein ID must exist in the protein FASTA.

### Protein FASTA

The first token of every protein FASTA header must equal the corresponding CDS `protein_id`:

```text
>NP_000001.1 protein description
MST...
```

Duplicate protein FASTA IDs cause a hard failure.

### CDS validity

`cubar::check_cds()` removes CDS records that are invalid under the configured genetic code. The original CDS remains represented in `gene_features.tsv`, but metrics are missing and `cds_qc_pass=false`.

---

## 3. Genetic codes

For genetic code `1`, tRNAscan-SE uses its standard code. For alternative codes, the workflow expects:

```text
resources/genetic_codes/trnascanse/<genetic_code>.gcode
```

The file should contain only codon reassignments required by tRNAscan-SE, for example:

```text
CTG Ser S
```

The tRNAscan-SE mode itself is configured globally by `trnascanse.domain` and defaults to `-E`; the sample-level `domain` column is used during tRNA-weight estimation.

## 4. KOfam resources

Required configuration:

```yaml
kofamscan:
  database:
    ko_list: resources/kofamscan/ko_list
    profiles_dir: resources/kofamscan/profiles
  ribosome_reference:
    brite_json: resources/kofamscan/ribosome_reference/ko03011.json
```

The profiles directory must contain one `Kxxxxx.hmm` file for every cytosolic eukaryotic ribosomal KO selected from the BRITE hierarchy.

## Optional external annotations

Set:

```yaml
paths:
  external_annotations: path/to/annotations.tsv
```

or leave it as `null`.

The table must contain `sample` and one unique join key:

- `protein_id`, or
- `seq_id`.

All additional columns are preserved. `go_terms` should contain GO IDs such as `GO:0006412`, separated by any text delimiter; valid IDs are extracted and normalized. Binary columns used by default statistics are:

```text
signal_peptide_present
tm_present
lcr_present
pfam_present
```

## 5. GO dictionary

The workflow accepts:

1. an existing valid TSV at `go_dictionary.path`, or
2. a local OBO file at `go_dictionary.obo_path`, or
3. a downloadable OBO URL at `go_dictionary.obo_url`.

The generated dictionary contains `go_terms`, `go_name` and `go_namespace`.

---

## 6. Example sample table

Example `config/samples_test.tsv`:

```tsv
sample	species	genetic_code	accession	domain	kingdom	phylum	lifestyle	genome	cds	proteome	include
Ztritici	Zymoseptoria_tritici	1	GCA_000219625.1	eukarya	Fungi	Ascomycota	plant_pathogen	GCA_000219625.1*.fna	GCA_000219625.1*.fna	GCA_000219625.1*.faa	True
Pplacenta	Postia_placenta	1	GCA_000006255.1	eukarya	Fungi	Basidiomycota	wood_saprotroph	GCA_000006255.1*.fna	GCA_000006255.1*.fna	GCA_000006255.1*.faa	True
Spombe	Schizosaccharomyces_pombe	1	GCA_000002945.2	eukarya	Fungi	Ascomycota	nectar_tap_saprotroph	GCA_000002945.2*.fna	GCA_000002945.2*.fna	GCA_000002945.2*.faa	True
```
---

## 7. Metadata used during aggregation and reports

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

## 8. Test data

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
