# Project map

## Main entry points

- `workflow/Snakefile` -- main Snakemake workflow.
- `config/config.yaml` -- production configuration.
- `config/config_test.yaml` -- test configuration.
- `config/samples.tsv` -- production sample metadata.
- `config/samples_test.tsv` -- test sample metadata.

## Rule modules

- `workflow/rules/trnascan_rule.smk` -- tRNA gene prediction.
- `workflow/rules/clean_tRNAscanSE_output.smk` -- raw tRNAscan-SE output conversion.
- `workflow/rules/prepare_trna_codon_counts_to_tai_rule.smk` -- amino acid–anticodon counting.
- `workflow/rules/codon_usage_metrics_rule.smk` -- R/cubar codon usage metrics.

## Scripts

- `workflow/scripts/convert_trnascanse_output_to_tsv.py` - conver unfriendly output to classic .tsv
- `workflow/scripts/prepare_trna_codon_counts_to_tai.py`
- `workflow/scripts/calculate_tAI.R`    - calculates tAI and other codon usage metrics
- `workflow/scripts/create_test_dataset.sh` -- downsampling

## Environments

- `workflow/envs/trnascanse.yaml`
- `workflow/envs/python.yaml`
- `workflow/envs/r.yaml`