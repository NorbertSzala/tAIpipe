# tRNAscan-SE genetic-code files

This directory contains genetic-code exception files for `tRNAscan-SE`.

`tRNAscan-SE` uses the option:

```bash
-g <file>
````

The file should contain codon reassignment exceptions relative to the standard genetic code.

## Format

```text
<Codon> <Three-letter amino acid> <One-letter amino acid>
```

Example:

```text
CTG     Ser     S
```

Use DNA codons with `T`, not RNA codons with `U`.

## Naming convention

Custom files are named by NCBI translation table ID:

```text
12.gcode
```

Copied upstream files are stored separately:

```text
copied_from_trnascanse/
```

## Standard code

For NCBI translation table 1, no `-g` file is used.

## Current custom files

| file       | NCBI table | meaning                                                           |
| ---------- | ---------: | ----------------------------------------------------------------- |
| `12.gcode` |         12 | Alternative Yeast Nuclear Code; `CTG` encodes Ser instead of Leu. |

