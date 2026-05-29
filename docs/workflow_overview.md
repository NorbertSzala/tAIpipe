
# Workflow overview

***

## 0. Input and files structure

### Structure
Below is the complete project directory map showing exactly where rule configurations, execution environments, scripts, data inputs, log traces, and per-sample execution outputs reside:

<!-- 
#TODO - to można dać jako osobny doc ... 
 -->

\```text
.
├── config/
│   ├── config.yaml                                 # Production workflow parameters
│   ├── config_test.yaml                            # Test profile configuration parameters
├── data/
│   ├── genetic_codes/                              # Global genetic code rules
│   │   ├── README.md
│   │   ├── genetic_codes.tsv
│   │   └── trnascanse/
│   │       ├── 12.gcode
│   │       ├── ...
│   └── tutorial_data/                              # Standardized tutorial workspace
│       ├── README.md
│       ├── input/
│       │   ├── CDS/                                # Input protein-coding sequence multi-FASTA files
│       │   │   ├── GCA_000002945.2_first_contig_CDS.fna
│       │   │   ├── ...
│       │   ├── genetic_codes.tsv
│       │   ├── genome/                             # Input truncated whole-genome assembly files
│       │   │   ├── GCA_000002945.2_first_contig_1Mbp.fna
│       │   │   ├── ...
│       │   ├── metadata/
│       │   │   ├── samples_test.tsv                # Active Snakemake compilation sample table
│       │   │   └── test_dataset.tsv                # Phenotypic descriptive metadata master sheet
│       │   └── trnascanse/                         # Local tRNAscan-SE genetic translation tables
│       │       ├── 12.gcode
│       │       ├── 3.gcode
│       │       ├── 4.gcode
│       │       ├── 6.gcode
│       │       ├── README.md
│       │       └── TRNASCANSE_VERSION.txt
│       └── output/                                 # Consolidated pipeline output directory
│           ├── logs/                               # Standard output (stdout) and error (stderr) traces
│           │   ├── Anidulans/
│           │   ├── ...
│           └── per_genome/                         # Localized, sample-isolated result tracks
│               ├── Anidulans/
│               │   ├── codon_metrics/              # R/cubar statistical outputs
│               │   │   ├── Anidulans_amino_acid_usage.csv
│               │   │   ├── Anidulans_cai.csv
│               │   │   ├── Anidulans_codon_counts.csv
│               │   │   ├── Anidulans_enc.csv
│               │   │   ├── Anidulans_fop.csv
│               │   │   ├── Anidulans_gc.csv
│               │   │   ├── Anidulans_gc3s.csv
│               │   │   ├── Anidulans_rscu.csv
│               │   │   ├── Anidulans_summary.tsv   # Consolidated index table per-sample
│               │   │   ├── Anidulans_tai.csv
│               │   │   └── Anidulans_trna_weights.csv
│               │   ├── counted_codons/             # Aggregated anticodon counts
│               │   └── trnascan/                   # Cleaned coordinate track and tRNA predictions
│               ├── Calbicans/                      # (Replicated structures for each sample name)
│               ├── ...
└── workflow/
    ├── Snakefile                                   # Master pipeline workflow engine
    ├── rules/                                      # Modular pipeline rules
    ├── schemas/                                    # Configuration structural validators
    ├── envs/                                       # Conda environment definition wrappers
    └── scripts/                                    # Internal script computational layer
\```

### Input

For an explicit explanation of input variables, column behaviors, and resolution rules, see the [Input Format Explanation](./input_format.md).

*** 

## 1. tRNAscan-SE

### Introduction & Purpose
The primary objective of this step is to map the entire transfer RNA (tRNA) repertoire of a given organism directly from its genomic assembly. tRNA tracking cannot be completed via blast-like sequence alignment alone because these non-coding RNAs rely heavily on highly conserved three-dimensional configurations (the classic cloverleaf structure) to interact with aminoacyl-tRNA synthetases and the ribosomally bound mRNA.

By scanning the entire genome sequence, the pipeline predicts:
* The chromosomal coordinates of every tRNA locus.
* The specific amino acid type it carries.
* The exact three-nucleotide **anticodon loop** sequence.
* Potential structural abnormalities, pseudo-elements, or variable intronic bounds.

The output generated here serves as the structural foundation for estimating cellular tRNA gene copy numbers, which act as a direct proxy for real-time cellular tRNA abundance.

### Execution Blueprint
* **Rule Name:** `run_trnascanse`
* **Rule Module Location:** `workflow/rules/trnascan_rule.smk`
* **Underlying Engine:** `tRNAscan-SE v2.0` (executed via the tool-specific environment `workflow/envs/trnascanse.yaml`).

### Input Profile
* **Sequence Input:** A resolved single-file genome assembly FASTA path extracted from `{DATA_GENOME}`.
* **Parameter Injections:** * `threads`: Extracted from `config.yaml` (`trnascanse -> threads`) to support scalable parallel processing.
  * `domain`: Mapped via the `samples.tsv` column. If `eukarya`, it passes `-E`; if `bacteria`, it passes `-B`; if `archaea`, it passes `-A`. This flag forces the tool to use distinct covariance structural baseline models optimized for that domain's specific tRNA features.

### Output Profile & Data Layout
`tRNAscan-SE` outputs a highly distributed collection of files written natively to `{PER_GENOME}/<sample>/trnascan/`:
1. `*_trnascan.out`: The core tabular textual report.
2. `*_trnascan.stats`: A high-level file containing alignment score cutoffs and a count summary of predicted tRNAs grouped by amino acid family.
3. `*_trnascan.bed` & `*_trnascan.gff`: Standard programmatic track formats used to view tRNA distributions inside genomic browsers (like IGV).
4. `*_trnascan.fasta`: A multi-FASTA file containing the raw predicted secondary-structure-validated sequences of the extracted tRNAs.

#### The Formatting Dilemma of Raw Output:
The primary output file (`*_trnascan.out`) uses a non-standard, legacy space-separated matrix decorated with multi-line text headers, trailing descriptive labels, and broken dashed separators (`--------`). This layout is fundamentally incompatible with automated data engineering tools like pandas or tidyverse. 

The raw structure looks like this:
\```text
Sequence   		tRNA    	Bounds  	tRNA	Anti	Intron Bounds	Inf
Name       	tRNA #	Begin   	End     	Type	Codon	Begin	End	Score	Note
--------   	------	-----   	------  	----	----─	----─	────	──────	──────
LR898468.1 	1	656159  	656232  	Ile	AAT	0	0	78.2
LR898468.1 	2	701452  	701523  	Arg	TCG	0	0	69.1   Pseudogene
\```

### Critical Biological Concept: The Anticodon vs Codon Mismatch
When reviewing the `Codon` or `Anti` columns inside `tRNAscan-SE` output, it is vital to note that **the software outputs DNA-triplet anticodons, not mRNA codons.**

1. **Transcription Alignment:** If the output lists `Anti` as `AAT`, it represents the genomic DNA coordinate of the tRNA gene. Transcribed into the functional tRNA molecule, this loop sequence becomes the RNA anticodon **`AAU`**.
2. **mRNA Target Pairing:** To understand which mRNA codon this specific tRNA molecule is built to target during translation, you must evaluate its **reverse complement**. For an RNA anticodon of `5'-AAU-3'`, the perfectly complementary mRNA codon sequence reads `3'-UUA-5'`, which matches the standard mRNA code codon **`ATT`** (or `AUU` in mature transcript text). 
3. This structural orientation translation is handled implicitly by the pipeline's downstream calculations to correctly match the tRNA pool with the coding sequence (CDS) tracks.


## 2. clean_tRNAscanSE_output

### Introduction & Purpose
Raw textual outputs generated by `tRNAscan-SE` contain multi-line headers, decorative alignment guides, and spaces that prevent direct parsing by data analysis tools. This process cleans and reformats the unaligned text matrix into a standardized, tabular data format. It evaluates records for minimum score requirements, filters biological features based on configuration variables, and establishes structured types for downstream routines.

### Execution Blueprint
* **Rule Name:** `clean_tRNAscanSE_output`
* **Module:** `workflow/rules/clean_tRNAscanSE_output.smk`
* **Script Component:** `workflow/scripts/convert_trnascanse_output_to_tsv.py`
* **Runtime Environment:** `workflow/envs/python.yaml`

### Operational Configuration Variables
The rule behavior is dictated by parameter adjustments inside `config/config.yaml`:
\```yaml
trnascanse_clean:
  keep_pseudo: false  # If false, removes all identified pseudo-tRNA genes from the output
  remove_nnn: true    # If true, discards entries with undetermined or unassigned anticodons
\```

### Target Routing
* **Input:** `{PER_GENOME}/<sample>/trnascan/<sample>_trnascan.out`
* **Output:** `{PER_GENOME}/<sample>/trnascan/<sample>_trnascan.tsv`

### Data Transformation Schema
The processing script applies structural validation rules:
1. Strips text headers and structural divider rows.
2. Formats all data rows against a structured list of output columns.
3. Asserts that the anticodon string matches a valid 3-letter IUPAC code configuration `^[ACGTN]{3}$`.
4. Asserts that structural validation alignment scores are positive numbers ($> 0.0$).
5. Flits logic flags to drop rows where `is_pseudo == True` or `anticodon == "NNN"` according to active settings.

The resulting clean tabular output uses the following standard configuration:
| Column Name | Description |
| :--- | :--- |
| `seq_name` | Host contig or chromosomal structural identifier string. |
| `trna_number` | Ordinal element tracking index relative to position within the host contig. |
| `begin` | Leftmost sequence coordinate index on the chromosome. |
| `end` | Rightmost sequence coordinate index on the chromosome. |
| `trna_type` | The assigned target amino acid passenger identity (e.g., `Leu`, `Ala`). |
| `anticodon` | The identified 3-nucleotide DNA-triplet anticodon sequence loop. |
| `score` | Calculated bit score evaluation from the structural covariance model. |
| `is_pseudo` | Boolean flag marking structural pseudo-elements lacking translation capability. |

***

## 3. prepare_trna_codon_counts_to_tai

### Introduction & Purpose
Downstream calculation engines do not require localization coordinate tracks or transcript sequence fragments. They require an explicit summary showing the total copy number concentration available for each unique anticodon. This module collapses the clean coordinate matrix into an aggregated count table to build an operational profile of the cell's transfer RNA gene pool.

### Execution Blueprint
* **Rule Name:** `prepare_trna_codon_counts_to_tai`
* **Module:** `workflow/rules/prepare_trna_codon_counts_to_tai_rule.smk`
* **Script Component:** `workflow/scripts/prepare_trna_codon_counts_to_tai.py`
* **Runtime Environment:** `workflow/envs/python.yaml`

### Target Routing
* **Input:** `{PER_GENOME}/<sample>/trnascan/<sample>_trnascan.tsv`
* **Output:** `{PER_GENOME}/<sample>/counted_codons/<sample>_aaa_counts.tsv`

### Aggregation Mechanism
1. The script reads the clean tabular data from the previous step.
2. It verifies that `trna_type` and `anticodon` are present in the table.
3. It filters out any remaining pseudo-tRNA elements where `is_pseudo == True`.
4. It concatenates the two tracking properties to create a unique identifier token (`anticodon_id`) matching the pattern: `[AminoAcid]-[Anticodon]` (e.g., `Ala-CGC`).
5. It runs a grouping routine that tallies total occurrences for each identifier.
6. The data is sorted alphabetically by `anticodon_id` and exported as a two-column distribution table.

#### Output Layout Example:
\```text
anticodon_id    count
Ala-CGC         4
Gly-TCC         7
Leu-CAG         3
\```


## 4. codon_usage_metrics

### Introduction & Purpose
This step serves as the core analytical engine of `tAIpipe`. It integrates the discrete tRNA availability profiles (gene copy number counts) generated in Step 3 with the actual codon demand derived from the coding sequences (CDS). 

By leveraging the computational layers of the R package `cubar`, this module maps tRNA abundance to codon-anticodon pairing rules (including domain-specific wobble interactions). It translates genomic abundance into relative adaptation weights for all 61 sense codons, projects these weights back onto every individual gene, and computes a comprehensive slate of evolutionary, compositional, and translational efficiency metrics.

### Execution Blueprint
* **Rule Name:** `codon_usage_metrics`
* **Rule Module Location:** `workflow/rules/codon_usage_metrics_rule.smk`
* **Script Component:** `workflow/scripts/calculate_tAI.R`
* **Underlying Engine:** Rscript running `cubar`, `Biostrings`, and `tidyverse` components (executed via the dedicated runtime environment `workflow/envs/r-cubar.yaml`).

### Input Profile
* **Tabular Input:** `{PER_GENOME}/<sample>/counted_codons/<sample>_aaa_counts.tsv` (The validated anticodon copy-number profile).
* **Sequence Input:** A resolved single-file protein-coding sequence multi-FASTA path from `{DATA_CDS}` (containing the nucleotide open reading frames of the target organism).
* **Parameter Injections:** * `genetic-code`: Extracted dynamically per sample from `samples_test.tsv` to ensure non-standard nuclear or mitochondrial translations are applied correctly.
  * `domain`: Extracted dynamically per sample from `samples_test.tsv` (`eukarya`, `bacteria`, or `archaea`). This parameter dictates the internal wobble-pairing coefficient matrix utilized during weights estimation.

### Output Profile & Data Layout

All results are systematically exported as tabular matrices into `{PER_GENOME}/<sample>/codon_metrics/`:



| File Name | Metric Name | Description |
| --- | --- | --- |
| `*_trna_weights.csv` | Relative adaptiveness weights ($w_i$) | The calculated relative adaptiveness weights for each of the 61 sense codons. |
| `*_tai.csv` | tRNA Adaptation Index | Gene-by-gene raw tRNA Adaptation Index scores. |
| `*_amino_acid_usage.csv` | Amino acid usage matrix | Global frequency distribution matrix of translated amino acids per gene. |
| `*_fop.csv` | Frequency of Optimal Codons | Frequency of Optimal Codons values calculated per gene sequence. |
| `*_cai.csv` | Codon Adaptation Index | Codon Adaptation Index scores mapped across all sequences. |
| `*_enc.csv` | Effective Number of Codons | Effective Number of Codons calculations evaluating local synonymous bias strength. |
| `*_rscu.csv` | Relative Synonymous Codon Usage | Relative Synonymous Codon Usage profiles for all synonymous groups. |
| `*_gc.csv` & `*_gc3s.csv` | Fractional GC and GC3s distribution | Fractional GC distribution vectors tracking total and third-position synonymous nucleotide bias. |
| `*_summary.tsv` | **Master consolidated summary table** | **The master single-sample consolidated summary table.** This file merges every single calculated metric into a unified rows-by-columns format where each row represents a unique gene ID from the input CDS file. |

***

### Detailed Metrics Reference
For an explicit explanation of metrics, their definitions and interpretations see the [Detailed Metrics Reference](./metrics.md).





