# Input format

This document provides reference for the core input files required to run `tAIpipe`. It details how files need to be structure and how define sample properties for statistical grouping.

***

## 0. Core Input Architecture

`tAIpipe` requires three types of inputs to compute codon usage and translational adaptation metrics:
1. **Genomic Sequence Files (`{DATA_GENOME}`)**
2. **Coding Sequence Files (`{DATA_CDS}`)**
3. **Tabular Control and Metadata Sheets (`metadata/dataset_test`)**

<!-- 
#TODO - Norbert - możesz dodać jak te tabel użyć :) 
 -->

### 1. Genome FASTA Files (`{DATA_GENOME}`)
* **What it is:** DNA sequences representing the whole-genome assemblies (chromosomes, scaffolds, or contigs) of the target organisms in standard FASTA format (`.fna` or `.fasta`).
* **Purpose:** tRNA genes are non-coding RNA (ncRNA) elements. They do not reside inside exons of protein-coding genes. Instead, they are distributed across intergenic regions and may contain unique intronic spaces. The complete structural context of the genome is mandatory for structural covariance scanning tools to detect them.

### 2. Coding Sequences (CDS) FASTA Files (`{DATA_CDS}`)
* **What it is:** Multi-FASTA files containing only the nucleotide sequences that directly code for proteins (from the initiation codon START to the termination codon STOP).
* **Purpose:** These sequences are parsed to count the exact frequencies of all 61 sense codons within the organism's expressed genes. These frequencies are later weighted against cellular tRNA abundance.


## 2. Tabular Management Files

The pipeline splits configuration into two files to separate **workflow orchestration** from **biological/ecological metadata grouping**.

### A. Experiment metadata - dataset (`data/tutorial_data/input/metadata/dataset_test.tsv`)
This is the active control file read directly by `workflow/Snakefile`. Its parameters determine exactly which jobs Snakemake compiles into the directed acyclic graph (DAG).

#### Required Schema Columns - dataset:

| Column Name | Data Type | Operational Category | Purpose & Workflow Impact |
| :--- | :--- | :--- | :--- |
| `sample` | String | Active Execution Parameter | Unique identifier string for the organism (e.g., `Spombe`). Used as the key wildcard value to instantiate isolated output file targets under `{PER_GENOME}/{sample}/`. |
| `species` | String | Downstream Analytical Grouping | The binomial nomenclature name of the organism (e.g., `Schizosaccharomyces_pombe`). <br> *Status:* `#TODO - Reserved for Downstream Aggregation and Final Report Formatting.` |
| `genetic_code` | Integer | Active Execution Parameter | The standard NCBI translation table numeric token (e.g., `1`, `3`, `12`). This variable is passed directly to computational parameters across software rules to determine active codon structures. |
| `domain` | String | Active Execution Parameter | Must strictly equal `eukarya`, `bacteria`, or `archaea`. Sets structural models for `tRNAscan-SE` (`-E`, `-B`, `-A`) and selects matching wobble pairing penalty matrices inside `cubar`. |
| `accession` | String | Downstream Analytical Grouping | The official NCBI assembly reference id (e.g., `GCA_000002945.2`). Used to cross-reference primary source archives. <br> *Status:* `#TODO - Reserved for Downstream Aggregation.` |
| `kingdom` | String | Downstream Analytical Grouping | Broad taxonomic kingdom classification (e.g., `Fungi`). Injected as a global filtering criterion during high-level analysis steps. <br> *Status:* `#TODO - Reserved for Downstream Aggregation.` |
| `phylum` | String | Downstream Analytical Grouping | Phylum-level taxonomy string used to cluster evolutionary features. <br> *Status:* `#TODO - Reserved for Downstream Aggregation.` |
| `lifestyle` | String | Downstream Analytical Grouping | Ecological niche category designation (e.g., `saprotroph`, `pathogen`). Used to evaluate selective adaptations in codon optimization. <br> *Status:* `#TODO - Reserved for Downstream Aggregation.` |
| `genome` | String | Active Execution Parameter | Filename token or glob/regex pattern pointing to the target sequence inside `{DATA_GENOME}` (e.g., `GCA_000002945.2*.fna`). |
| `cds` | String | Active Execution Parameter | Filename token or glob/regex pattern pointing to the target protein-coding library inside `{DATA_CDS}`. |
| `bed` | String | Developmental / Test Artifact | Pointer to an optional coordinate mapping file (`.bed`) defining specific target regions. <br> *Status:* `#TODO - Reserved for Downstream Aggregation.` |
| `downsampled_fasta` | String | Developmental / Test Artifact | Reference name for custom fast-testing inputs containing pre-selected genomic segments. <br> *Status:* `#TODO - Reserved for Downstream Aggregation.` |
| `first_contig` | String | Developmental / Test Artifact | Tracking label indicating the target file holding the isolated first contig sequence of the organism assembly. <br> *Status:* `#TODO - Reserved for Downstream Aggregation.` |
| `include` | Boolean | Active Execution Parameter | Logic gate evaluating to `True` or `False`. Rows matching `False` are culled during initialization, entirely excluding the target sample from processing. |



### B. Report metadata -  (`data/tutorial_data/input/metadata/samples_test.tsv`)
#TODO - to dopiero zrobimy ale masz wzorzec jak to będzie wyglądać 
#TODO - dodać że się to tworzy samo
* **What it is:** A comprehensive matrix describing the biological taxonomy (`kingdom`, `phylum`, `class`) and ecological characteristics (`lifestyle`, `microenvironment`) of every available sample.
* **Workflow Impact:** This file is completely ignored during the initial heavy computational rules (`run_trnascanse`, `codon_usage_metrics`). It is injected exclusively during **Step 5 (Aggregation and Report)**, allowing the pipeline to automatically group computed tAI, FOP, and ENC metrics by phenotypes—such as comparing whether wood-decaying saprotrophs exhibit higher translational efficiency for certain gene sets compared to plant pathogens.

***


### Matching Constraints:
1. **Zero Files Found:** If a pattern (e.g., `GCA_99999*.fna`) fails to match any file inside `{DATA_GENOME}`, the engine triggers a `MissingInputException`.
2. **Multiple Files Found:** If a pattern matches more than one file (e.g., `GCA_0002*.fna`) , the pipeline intentionally crashes instantly with a `Multiple files found for pattern` error. **Each pattern must resolve to exactly one unique sequence file.**

***
<!-- 
DEPRACATED - nie posiadamy skryptu który działa na dowolnych danach / ścieżce 
## 4. Test Data Generation and Downsampling

For prototyping and testing local execution configurations, the workflow includes an automated downsampling bash script.

* **Executing Script:** `workflow/scripts/create_test_dataset.sh`
* **Mechanics:** 1. The script isolates the very first contig or chromosome from a heavy production genome FASTA file.
  2. It truncates that contig precisely at the first 1,000,000 bp (controlled via `-L` parameter).
  3. It searches the matching master CDS file and extracts up to 1,000 records that map back to that specific contig ID, outputting a fully functional, lightweight mock environment inside `resources/test_data/`.
 -->
