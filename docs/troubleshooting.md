# Troubleshooting

## `MissingInputException`

Check whether file patterns in `config/samples.tsv` resolve to existing files.

```bash
ls data/genome/
ls data/CDS/
```

Each pattern should match exactly one file.

## `Multiple files found for pattern`

The workflow intentionally fails if a glob pattern matches more than one file. Replace broad patterns such as:

**GCA_000002945.2*.fna**

with a more specific filename.

## `ModuleNotFoundError: pandas`

Make sure the workflow is executed with the correct conda environment:

```bash
snakemake --profile workflow/profiles/default --use-conda
```

or install manually:

```bash
micromamba install -n taipipe -c conda-forge pandas
```

## tRNAscan-SE not found

Use the provided conda environment:

```bash
snakemake --profile workflow/profiles/default --use-conda
```

or check installation:

```bash
tRNAscan-SE --version
```

## Snakemake lock error

Unlock the working directory:

```bash
snakemake --unlock
```
