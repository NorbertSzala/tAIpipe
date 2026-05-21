#!/usr/bin/env bash

set -euo pipefail

snakemake -s workflow/Snakefile \
    --configfile config/config_test.yaml \
    --cores 2 \
    --use-conda \
    --dry-run