# Output format

`tAIpipe` writes results per genome and aggregated outputs.

## Per-genome output

Default location:

```text
results/per_genome/<sample>/
```

Example:

```txt
results/per_genome/Spombe/
├── trnascan/
│   ├── Spombe_trnascan.out
│   ├── Spombe_trnascan.tsv
│   ├── Spombe_trnascan.stats
│   ├── Spombe_trnascan.bed
│   ├── Spombe_trnascan.gff
│   └── Spombe_trnascan.fasta
└── counted_codons/
    └── Spombe_aaa_counts.tsv
```

## tRNA anticodon count table

Produced by:

```bash
workflow/scripts/prepare_trna_codon_counts_to_tai.py
```

Format:

| column              | description                                                 |
| ------------------- | ----------------------------------------------------------- |
| anticodon_id	| amino acid and anticodon pair, e.g. Ala-CGC   |
| count |	number of predicted tRNA genes with this amino acid–anticodon pair  |


Example:

| anticodon_id              | count                                                 |
| ------------------- | ----------------------------------------------------------- |
| Ala-CGC	| 4 |
| Gly-TCC	| 7 |
| Leu-CAG	| 3 |


## Aggregated results

Default location:

```txt
results/aggregated/
├── tables/
├── plots/
└── reports/
```