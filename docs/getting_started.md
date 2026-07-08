# Getting started

## 1. Prepare the repository

Run all commands from the repository root.

```bash
cd tAIpipe
```

The recommended backend is Conda through a Snakemake profile. Environments are defined in `workflow/envs/`.

## 2. Provide required resources

Production execution requires:

```text
<CDS directory>/
<genome directory>/
<proteome directory>/
resources/genetic_codes/trnascanse/
resources/kofamscan/ko_list
resources/kofamscan/profiles/Kxxxxx.hmm
resources/kofamscan/ribosome_reference/ko03011.json
```

Set the first three directories in `config/config.yaml`:

```yaml
paths:
  data_cds: /path/to/CDS
  data_genome: /path/to/genomes
  data_proteome: /path/to/proteomes
```

The GO OBO file is downloaded automatically only when the configured local file and output dictionary are missing.

## 3. Configure samples

Edit `config/samples.tsv`. Every row with `include=true` becomes one workflow sample. Each `genome`, `cds` and `proteome` pattern must resolve to exactly one file in its configured directory.

Minimal example:

```tsv
sample	species	genetic_code	accession	domain	kingdom	phylum	lifestyle	genome	cds	proteome	include
Spombe	Schizosaccharomyces_pombe	1	GCA_000002945.2	eukarya	Fungi	Ascomycota	saprotroph	GCA_000002945.2*.fna	GCA_000002945.2*.fna	GCA_000002945.2*.faa	True
```

See [`input_format.md`](input_format.md) before adding a new dataset.

## 4. Validate configuration and DAG

Test profile:

```bash
snakemake -n --profile workflow/profiles/test
```

Production profile:

```bash
snakemake -n --profile workflow/profiles/production
```

A successful dry run confirms that configuration validation, sample selection and file dependencies can be resolved. It does not validate biological content or execute tools.

## 5. Run the workflow

```bash
snakemake --profile workflow/profiles/test
```

or:

```bash
snakemake --profile workflow/profiles/production
```

Useful flags for manual runs:

```bash
--rerun-incomplete
--printshellcmds
--show-failed-logs
--keep-going
```

## 6. Check completion

The workflow is complete when all targets from `rule all` exist:

```text
results/tables/gene_features.tsv
results/tables/genome_summary.tsv
results/tables/codon_profiles.tsv
results/statistics/gene_feature_tests.tsv
results/statistics/genome_group_tests.tsv
results/statistics/go_enrichment.tsv
```

Per-sample logs are stored under `logs/<sample>/`; aggregate steps use `logs/tables/`, `logs/statistics/` or `logs/aggregated/`.

## 7. Re-run selected stages

Examples:

```bash
snakemake --profile workflow/profiles/test build_gene_features
snakemake --profile workflow/profiles/test build_genome_summary
snakemake --profile workflow/profiles/test compute_statistics
```

Force one target only when its inputs or implementation changed:

```bash
snakemake --profile workflow/profiles/test -R build_gene_features
```
