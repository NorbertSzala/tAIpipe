# Test-data support and current test status

This repository currently provides a deterministic test-data generator and a
Snakemake test profile. It does not yet contain implemented unit or integration
test cases: `tests/unit/` and `tests/integration/` are empty placeholders.

The generator is useful for producing a small, internally linked genome/CDS/
protein/GFF fixture from real assemblies. It is not itself evidence that the
complete workflow passes.

## What the generator creates

`workflow/scripts/preprocessing/create_test_dataset.sh` processes each selected
accession independently:

1. It requires exactly one matching genome FASTA, CDS FASTA, protein FASTA and
   GFF/GFF3 file in the supplied source directories.
2. It retains the first genome contig, truncated to `-L` bases (default:
   1,000,000).
3. It selects CDS records assigned to that contig and contained within the
   retained interval. If none can be matched, it falls back to the first `-N`
   CDS records (default: 1,000) and reports that fallback.
4. It extracts `[protein_id=...]` values from the selected CDS headers and
   requires every selected protein ID to exist in the protein FASTA.
5. It writes only matching protein records and GFF features located within the
   retained interval.

The output is:

```text
tests/data/
|-- genome/
|-- CDS/
|-- proteome/
`-- gff/
```

The `-O` option changes the output root. The generated genome filename retains
the `_first_contig_1Mbp` suffix even when a non-default `-L` value is used; the
actual retained length is printed by the script.

## Required invocation

Run the Bash script from the repository root and provide all four source data
directories:

```bash
bash workflow/scripts/preprocessing/create_test_dataset.sh \
  -C /path/to/CDS \
  -G /path/to/genomes \
  -P /path/to/proteomes \
  -F /path/to/gff \
  -A "GCA_000219625.1,GCA_000006255.1" \
  -O tests/data \
  -L 1000000 \
  -N 1000
```

Arguments:

| Option | Meaning |
| --- | --- |
| `-C` | Source CDS FASTA directory |
| `-G` | Source genome FASTA directory |
| `-P` | Source protein FASTA directory |
| `-F` | Source GFF/GFF3 directory |
| `-A` | Comma- or space-separated assembly accessions |
| `-O` | Output root; default `tests/data` |
| `-L` | Maximum number of bases retained from the first contig |
| `-N` | Maximum fallback CDS count when contig matching fails |

Compressed `.gz` inputs are accepted. The script rejects zero or multiple
matches for an accession because silently choosing one file could break the
CDS-protein mapping.

## Connecting fixtures to Snakemake

`workflow/profiles/test/config.yaml` points to `config/config_test.yaml`, which
in turn points to `tests/data/` and `config/samples_test.tsv`. After generating
fixtures, the included sample rows and filename patterns must match the files
under `tests/data/` exactly.

Current limitation: `config/config_test.yaml` is not synchronized with the
current `workflow/schemas/config.schema.yaml`. It lacks the required top-level
sections `gene_protein_map`, `plots`, `script_suggestion_plots`,
`codon_profile_plots`, `go_plots`, `pfam_lcr_plots` and `plot_style`; several
existing sections also lack required nested fields. Therefore:

```bash
snakemake --dry-run --profile workflow/profiles/test
```

is expected to fail configuration validation until the test config is updated.
The test profile must not be presented as a working integration test or release
gate in its current state.

Once the test configuration is synchronized, the minimum meaningful checks are
a full dry-run followed by a run against `tests/data/`, with inspection of the
tRNA, metric and KOfam QC tables. Automated assertions for those outputs still
need to be implemented under `tests/unit/` and `tests/integration/`.
