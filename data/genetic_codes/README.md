# Genetic code resources

This directory stores genetic-code resources used by `tAIpipe`.

## Directory layout

```text
data/genetic_codes/
├── genetic_codes.tsv
├── tables/
└── trnascanse/
```

## Files

| path                | purpose                                                                         |
| ------------------- | ------------------------------------------------------------------------------- |
| `genetic_codes.tsv` | Manifest mapping NCBI translation table IDs to local genetic-code resources.    |
| `tables/`           | Full codon translation tables for CDS/codon-usage analyses. Currently optional. |
| `trnascanse/`       | Genetic-code exception files compatible with `tRNAscan-SE -g`.                  |

## tRNAscan-SE genetic-code files

`tRNAscan-SE` does not accept NCBI genetic-code numbers directly.

The `-g` option expects a file containing codon reassignment exceptions relative to the standard code.

Format:

```text
<Codon> <Three-letter amino acid> <One-letter amino acid>
```

Example for NCBI translation table 12, Alternative Yeast Nuclear Code:

```text
CTG     Ser     S
```

For standard nuclear code, no `-g` argument should be passed.

## Current policy

For nuclear fungal genomes in this project:

| organism type                         | `genetic_code` in `samples.tsv` | tRNAscan-SE `-g`                         |
| ------------------------------------- | ------------------------------: | ---------------------------------------- |
| Standard-code fungi                   |                               1 | not used (default)                       |
| `Candida albicans` / CTG-clade yeasts |                              12 | `data/genetic_codes/trnascanse/12.gcode` |
