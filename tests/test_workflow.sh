#!/usr/bin/env bash

set -euo pipefail

snakemake -n --profile workflow/profiles/test
snakemake --profile workflow/profiles/test --rerun-incomplete --printshellcmds