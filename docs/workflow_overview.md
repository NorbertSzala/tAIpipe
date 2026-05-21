# Workflow overview

`tAIpipe` is a Snakemake workflow for comparative tRNA adaptation and codon usage analysis.

## Main steps

1. Read sample metadata from `config/samples.tsv`
2. Resolve genome and CDS input files
3. Run `tRNAscan-SE` on genome or downsampled genome FASTA
4. Convert raw `tRNAscan-SE` output into clean TSV
5. Count amino acid–anticodon pairs
6. Calculate tAI and codon usage metrics
7. Export per-genome and aggregated result tables

## Core rules

| rule | purpose |
|---|---|
| `run_trnascanse` | predicts tRNA genes |
| `clean_tRNAscanSE_output` | converts raw tRNAscan-SE output to TSV |
| `prepare_trna_codon_counts_to_tai` | counts amino acid–anticodon pairs |
| `codon_usage_metrics` | calculates tAI, CAI, ENC, RSCU, GC, GC3s |

## Minimal run

```bash
snakemake --profile workflow/profiles/default
```

## Dry run

```bash
snakemake --profile workflow/profiles/default --dry-run
```