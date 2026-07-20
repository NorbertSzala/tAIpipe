# tRNAscan-SE genetic-code exceptions

This directory contains the genetic-code exception files passed to
`tRNAscan-SE` with `-g`. They change tRNA isotype assignment for translation
tables that differ from the standard code.

## How the workflow selects a file

For every included sample, `workflow/common/03_input_functions.smk` reads the
integer `genetic_code` from the sample sheet.

- Code `1` uses the tRNAscan-SE default and no `-g` argument.
- Any other code requires
  `resources/genetic_codes/trnascanse/<genetic_code>.gcode`.
- A missing file stops DAG construction instead of silently applying the
  standard code.

The same sample-level genetic-code ID is later used by the R metric and table
scripts. Keeping both interpretations synchronized is necessary because code
changes affect tRNA isotypes, stop codons, synonymous families, RSCU, CAI and
tAI.

The tRNAscan-SE domain mode is independent of this file. It is selected from
the sample's `domain` field (`Eukarya`, `Bacteria` or `Archaea`).

## File format

Each non-comment line contains a DNA codon or accepted IUPAC codon pattern, a
three-letter amino-acid name and a one-letter amino-acid code:

```text
<DNA codon/pattern> <three-letter amino acid> <one-letter amino acid>
```

Use `T`, not `U`. The files specify only reassignments relative to the standard
code; they are not complete codon tables.

Example:

```text
CTG Ser S
```

## Files currently provided

| File | NCBI translation table | Exceptions encoded in the file |
| --- | ---: | --- |
| `3.gcode` | 3, Yeast Mitochondrial | `TGA -> Trp`, `ATA -> Met`, `CTN -> Thr` |
| `4.gcode` | 4, Mold/Protozoan/Coelenterate Mitochondrial and Mycoplasma/Spiroplasma | `TGA -> Trp` |
| `6.gcode` | 6, Ciliate/Dasycladacean/Hexamita Nuclear | `TAR -> Gln` |
| `12.gcode` | 12, Alternative Yeast Nuclear | `CTG -> Ser` |

The sample schema accepts additional NCBI translation-table IDs, but they are
not usable with tRNAscan-SE until a matching reviewed `.gcode` file is added to
this directory. Adding a filename without validating the actual reassignment is
not sufficient.
