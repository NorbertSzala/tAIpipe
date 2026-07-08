# Output files

Output roots are controlled by `paths` in the selected config file. The examples below use production defaults.

## Directory layout

```text
results/
├── per_genome/<sample>/
│   ├── trnascan/
│   ├── counted_codons/
│   ├── qc/
│   ├── tables/
│   ├── kofamscan/
│   └── codon_metrics/
├── tables/
│   ├── gene_features.tsv
│   ├── genome_summary.tsv
│   └── codon_profiles.tsv
└── statistics/
    ├── gene_feature_tests.tsv
    ├── genome_group_tests.tsv
    └── go_enrichment.tsv
```

Logs and benchmarks are written separately to `logs/` and `benchmarks/`.

## Per-genome outputs

### tRNAscan-SE

| File                      | Description                          |
| ------------------------- | ------------------------------------ |
| `<sample>_trnascan.out`   | Main raw tRNAscan-SE report          |
| `<sample>_trnascan.stats` | tRNAscan-SE summary statistics       |
| `<sample>_trnascan.bed`   | Genomic intervals                    |
| `<sample>_trnascan.gff`   | GFF annotations                      |
| `<sample>_trnascan.fasta` | Predicted tRNA sequences             |
| `<sample>_trnascan.tsv`   | Parsed and filtered prediction table |

### tRNA profile and QC

`counted_codons/<sample>_aaa_counts.tsv`:

| Column         | Meaning                                             |
| -------------- | --------------------------------------------------- |
| `anticodon_id` | Amino-acid and anticodon key, for example `Ala-CGC` |
| `count`        | Number of retained tRNA genes                       |

`qc/<sample>_trna_profile_qc.tsv` contains one row. Important columns:

| Column                | Meaning                                            |
| --------------------- | -------------------------------------------------- |
| `n_total_trnas`       | Sum of all valid counts                            |
| `n_elongator_trnas`   | Sum after excluding initiator/undetermined records |
| `n_unique_anticodons` | Number of distinct elongator anticodons            |
| `n_amino_acids`       | Number of represented standard amino acids         |
| `qc_status`           | `PASS`, `WARN` or `FAIL`                           |
| `qc_reasons`          | Semicolon-separated reasons                        |

### Identifier mapping

`tables/gene_protein_map.tsv`:

```text
sample, gene_id, gene_id_source, locus_tag, gene_name,
protein_id, cds_id, present_in_proteome, cds_description
```

### KofamScan ribosomal reference

| File                               | Main columns or content                                                                         |
| ---------------------------------- | ----------------------------------------------------------------------------------------------- |
| `ribosome_detail.tsv`              | Raw KofamScan `detail-tsv` output                                                               |
| `ribosome_significant_hits.tsv`    | `sample`, `gene_id`, `cds_id`, `protein_id`, `ko`, `threshold`, `score`, `evalue`, `definition` |
| `ribosome_gene_annotations.tsv`    | Gene-level KO assignments and ambiguity flags                                                   |
| `ribosomal_reference_gene_ids.txt` | One selected gene ID per line                                                                   |
| `ribosomal_reference_cds_ids.txt`  | One exact selected CDS ID per line                                                              |
| `ribosomal_reference_cds.fna`      | Nucleotide reference CDS FASTA used for CAI                                                     |
| `ribosome_qc.tsv`                  | KO coverage, hit counts and mapping summary                                                     |

Important `ribosome_qc.tsv` fields include `n_reference_kos`, `n_observed_kos`, `ko_coverage`, `n_ribosomal_genes`, `n_reference_cds` and `n_multi_ko_proteins`.

### Codon metrics

| File                            | Shape                 | Content                                                        |
| ------------------------------- | --------------------- | -------------------------------------------------------------- |
| `<sample>_codon_counts.csv`     | gene × codon          | Codon count matrix                                             |
| `<sample>_enc.csv`              | gene vector           | `seq_id`, `value`                                              |
| `<sample>_rscu.csv`             | codon table           | Genome-wide RSCU and CAI weight fields                         |
| `<sample>_reference_rscu.csv`   | codon table           | Ribosomal-reference RSCU and CAI weights                       |
| `<sample>_cai.csv`              | gene vector           | `seq_id`, `value`                                              |
| `<sample>_trna_weights.csv`     | codon table           | At least `codon`, `anticodon`, `trna_id`, `ac_level`, `W`, `w` |
| `<sample>_tai.csv`              | gene vector           | `seq_id`, `value`                                              |
| `<sample>_amino_acid_usage.csv` | gene × amino acid     | Amino-acid usage matrix                                        |
| `<sample>_fop.csv`              | gene vector           | `seq_id`, `value`                                              |
| `<sample>_gc.csv`               | gene vector           | `seq_id`, `value`                                              |
| `<sample>_gc3s.csv`             | gene vector           | `seq_id`, `value`                                              |
| `<sample>_summary.tsv`          | one row per valid CDS | `seq_id`, `ENC`, `CAI`, `FOP`, `GC`, `GC3s`, `tAI`             |

Package-generated codon tables may contain additional annotation columns such as `aa_code`, `amino_acid` and `subfam`.

## Canonical tables

### `gene_features.tsv`

**Key:** `sample + seq_id`  
**Unit:** one original CDS record.

Core column groups:

| Group               | Columns                                                                                                            |
| ------------------- | ------------------------------------------------------------------------------------------------------------------ |
| Sample metadata     | `sample`, `species`, `accession`, `domain`, `kingdom`, `phylum`, `lifestyle`, `genetic_code`                       |
| Identifiers         | `seq_id`, `gene_id`, `protein_id`, `uniprot_id`, `protein_name`                                                    |
| Sequence properties | `cds_length_nt`, `protein_length_aa`, `log_protein_length_aa`, `cds_length_multiple_of_three`, `has_terminal_stop` |
| QC                  | `cds_qc_pass`, `metrics_available`, `trna_qc_status`, `trna_qc_reasons`, `trna_qc_pass`                            |
| tRNA summary        | `n_total_trnas`, `n_elongator_trnas`, `n_unique_anticodons`, `n_trna_amino_acids`                                  |
| Metrics             | `ENC`, `ENC_expected`, `delta_ENC`, `CAI`, `FOP`, `tAI`, `tAI_z`, `tAI_percentile`, `GC`, `GC3s`                   |
| Default annotations | `signal_peptide_present`, `tm_present`, `lcr_present`, `pfam_present`, `pfam_terms`, `go_terms`                    |

Additional columns from the optional annotation table are preserved.

### `genome_summary.tsv`

**Key:** `sample`  
**Unit:** one included genome.

Main fields:

- sample metadata;
- `genome_file`, `n_contigs`, `genome_size_bp`, `genome_gc`, `genome_n_fraction`;
- `n_genes`, `n_genes_with_metrics`, metric/QC fractions;
- means and medians of tAI, ENC, CAI, FOP, GC, GC3s and `delta_ENC`;
- fractions of annotated gene features and GO coverage;
- tRNA QC fields prefixed with `trna_`;
- KofamScan QC fields prefixed with `kofam_`.

### `codon_profiles.tsv`

**Key:** `sample + codon`  
**Unit:** one DNA codon in one genome.

```text
sample, species, accession, domain, kingdom, phylum, lifestyle,
genetic_code, codon, is_stop_codon, aa_code, amino_acid, subfam,
codon_count, codon_frequency, genome_RSCU, reference_RSCU,
reference_CAI_weight, reference_optimal_codon, trna_anticodon,
trna_id, trna_copy_number, trna_absolute_weight, tRNA_weight
```

Stop codons are retained and marked by `is_stop_codon`; synonymous/tRNA fields may be missing for them.

## Statistical outputs

### `gene_feature_tests.tsv`

One row per configured binary feature:

```text
analysis, feature, term, estimate, std_error, statistic, p_value,
conf_low, conf_high, n_genes, n_genomes, covariates,
model_formula, status, q_value
```

### `genome_group_tests.tsv`

One row per configured metric and group variable:

```text
analysis, metric, group_variable, test, comparison, effect,
statistic, p_value, n_genomes, n_groups, status, q_value
```

### `go_enrichment.tsv`

One row per tested GO term and tAI tail:

```text
tail, go_id, go_name, go_namespace, n_genomes_with_term,
n_informative_genomes, total_genes_with_term, n_genomes_enriched,
n_genomes_depleted, n_genomes_neutral, common_odds_ratio,
conf_low, conf_high, statistic, p_value, status, q_value
```

An empty table with headers is valid when no GO term passes the configured support filters.
