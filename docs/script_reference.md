# Active script reference

This document covers every custom script used by the active numbered rules. External programs are listed separately.

## Resource preparation

### `workflow/scripts/resources/parse_go_obo_to_tsv.py`

**Rule:** `prepare_go_dictionary`

**Input**

- GO OBO file from `--input-obo`; downloaded from `--obo-url` only if missing.
- Existing output TSV may be reused when it already contains the required columns.

**Output**

- GO dictionary TSV with columns: `go_terms`, `go_name`, `go_namespace`.

**Purpose**

Removes obsolete terms and creates the compact lookup used to label GO enrichment results.

## tRNA preprocessing and QC

### `workflow/scripts/preprocessing/convert_trnascanse_output_to_tsv.py`

**Rule:** `convert_trnascanse_output_to_tsv`

**Input**

- Raw tRNAscan-SE `.out` file.
- Flags controlling pseudogene retention and removal of `NNN` anticodons.

**Output**

- Clean TSV with: `seq_name`, `trna_number`, `begin`, `end`, `trna_type`, `anticodon`, `intron_begin`, `intron_end`, `score`, `note`, `strand`, `start`, `stop`, `has_intron`, `is_pseudo`.

**Purpose**

Parses the fixed-width-like tRNAscan-SE report, standardizes coordinates and creates explicit strand/QC fields.

### `workflow/scripts/preprocessing/prepare_trna_codon_counts_to_tai.py`

**Rule:** `prepare_trna_codon_counts_to_tai`

**Input**

- Clean tRNAscan-SE TSV containing at least `trna_type` and `anticodon`.

**Output**

- Two-column TSV: `anticodon_id`, `count`.
- `anticodon_id` format: `<amino-acid>-<DNA anticodon>`, for example `Ala-CGC`.

**Purpose**

Removes pseudogenes, initiator/undetermined records and invalid anticodons, then counts elongator tRNA gene copies.

### `workflow/scripts/qc/validate_trna_profile.py`

**Rule:** `validate_trna_profile`

**Input**

- Two-column tRNA copy-number table.
- Configured minimum total tRNAs, unique anticodons and represented amino acids.

**Output**

- One-row QC TSV containing counts, structural-error counters, thresholds, `qc_status` and `qc_reasons`.

**Purpose**

Detects malformed identifiers, invalid counts, duplicates and biologically sparse profiles before tAI calculation.

## Identifier mapping and ribosomal reference

### `workflow/scripts/kofamscan/build_gene_protein_map.py`

**Rule:** `build_gene_protein_map`

**Input**

- NCBI-style CDS FASTA.
- Protein FASTA.
- Sample ID and strict matching flag from Snakemake.

**Output**

- `gene_protein_map.tsv` with: `sample`, `gene_id`, `gene_id_source`, `locus_tag`, `gene_name`, `protein_id`, `cds_id`, `present_in_proteome`, `cds_description`.

**Purpose**

Creates the explicit ID bridge required to convert protein-level KofamScan hits to gene and nucleotide CDS identifiers.

### `workflow/scripts/kofamscan/prepare_ribosomal_reference.py`

**Rule:** `prepare_kofam_ribosomal_reference`

**Input**

- KEGG BRITE `ko03011` JSON.
- KOfam `ko_list`.
- KOfam profile directory.

**Output**

- Ribosomal KO table with: `ko`, `subunit`, `brite_entry`, `threshold`, `threshold_available`, `score_type`, `profile_type`, `definition`.
- HAL profile-list file containing relative paths to selected `Kxxxxx.hmm` files.

**Purpose**

Restricts KofamScan to cytosolic eukaryotic large- and small-subunit ribosomal proteins and verifies profile availability.

### `workflow/scripts/kofamscan/parse_kofamscan_ribosome.py`

**Rule:** `parse_kofamscan_ribosome`

**Input**

- KofamScan detail TSV.
- Protein FASTA used as the query.
- Gene–protein–CDS mapping table.
- Allowed ribosomal KO table.

**Output**

- `ribosome_significant_hits.tsv`: one significant protein–KO hit.
- `ribosome_gene_annotations.tsv`: annotations collapsed to gene level.
- `ribosomal_reference_gene_ids.txt`: one gene ID per line.
- `ribosomal_reference_cds_ids.txt`: one exact CDS FASTA ID per line.
- `ribosome_qc.tsv`: coverage and mapping summary.

**Purpose**

Retains significant restricted KO assignments, validates query IDs and maps every hit to gene and CDS identifiers.

### `workflow/scripts/kofamscan/extract_reference_cds.py`

**Rule:** `extract_ribosomal_reference_cds`

**Input**

- Full nucleotide CDS FASTA.
- Exact CDS identifier list from the KofamScan parser.

**Output**

- FASTA containing only selected cytosolic ribosomal CDS records.

**Purpose**

Creates the sample-specific reference set used to estimate CAI weights.

## Metric calculation

### `workflow/scripts/metrics/calculate_tAI.R`

**Rule:** `codon_usage_metrics`

**Input**

- Full CDS FASTA.
- Ribosomal reference CDS FASTA.
- tRNA copy-number TSV.
- Sample ID, NCBI genetic code, biological domain and minimum reference size.

**Output**

- `*_codon_counts.csv`
- `*_enc.csv`
- `*_rscu.csv`
- `*_reference_rscu.csv`
- `*_cai.csv`
- `*_trna_weights.csv`
- `*_tai.csv`
- `*_amino_acid_usage.csv`
- `*_fop.csv`
- `*_gc.csv`
- `*_gc3s.csv`
- `*_summary.tsv`

**Purpose**

Filters invalid CDS using the selected genetic code, estimates reference-derived CAI weights, estimates domain-specific tRNA weights and calculates all gene-level codon metrics through `cubar`.

## Canonical table construction

### `workflow/scripts/tables/build_gene_features.R`

**Rule:** `build_gene_features`

**Input**

- Sample table.
- Per-sample `*_summary.tsv` metric files.
- Original CDS FASTA files.
- Per-sample tRNA QC tables.
- Optional external annotation TSV.

**Output**

- `gene_features.tsv`, one row per original CDS record.

**Purpose**

Combines sequence/header metadata, codon metrics, taxonomy, QC and annotations; derives `log_protein_length_aa`, expected ENC, `delta_ENC`, within-genome `tAI_z` and `tAI_percentile`.

### `workflow/scripts/tables/build_genome_summary.R`

**Rule:** `build_genome_summary`

**Input**

- `gene_features.tsv`.
- Sample table.
- Genome FASTA files.
- Per-sample tRNA and KofamScan QC files.

**Output**

- `genome_summary.tsv`, exactly one row per included sample.

**Purpose**

Aggregates gene metrics and annotation fractions, calculates genome size/GC/N content and appends QC summaries.

### `workflow/scripts/tables/build_codon_profiles.R`

**Rule:** `build_codon_profiles`

**Input**

- Per-sample codon-count, genome RSCU, reference RSCU and tRNA-weight files.
- Sample table.

**Output**

- `codon_profiles.tsv`, one row per sample and DNA codon.

**Purpose**

Creates a comparable genome-by-codon table containing codon demand, synonymous usage, CAI reference weights and tRNA adaptation weights.

## Statistics

### `workflow/scripts/statistics/compute_statistics.R`

**Rule:** `compute_statistics`

**Input**

- `gene_features.tsv`.
- `genome_summary.tsv`.
- Configured binary features, gene covariates, genome metrics, group variables and FDR method.

**Output**

- `gene_feature_tests.tsv`: one mixed-model result per binary feature.
- `genome_group_tests.tsv`: one non-parametric test per metric and grouping variable.

**Purpose**

Tests gene-level feature associations with within-genome standardized tAI and genome-level differences between taxonomy/lifestyle groups.

### `workflow/scripts/statistics/compute_go_enrichment_cmh.R`

**Rule:** `compute_go_enrichment_cmh`

**Input**

- `gene_features.tsv` containing `sample`, `gene_id`, finite `tAI` and `go_terms`.
- GO dictionary.
- Tail-size, minimum-support and FDR parameters.

**Output**

- `go_enrichment.tsv`, one row per tested GO term and high/low-tAI tail.

**Purpose**

Selects high- and low-tAI genes independently within each genome and performs a Cochran–Mantel–Haenszel enrichment test stratified by genome.

## External executables

| Program           | Rule                     | Input                                    | Output                                                 |
| ----------------- | ------------------------ | ---------------------------------------- | ------------------------------------------------------ |
| `tRNAscan-SE`     | `run_trnascanse`         | Genome FASTA                             | `.out`, `.stats`, `.bed`, `.gff`, predicted tRNA FASTA |
| `exec_annotation` | `run_kofamscan_ribosome` | Protein FASTA, restricted HAL, `ko_list` | KofamScan detail TSV                                   |
