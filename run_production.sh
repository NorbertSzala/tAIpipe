#!/usr/bin/env bash
set -euo pipefail

snakemake --profile workflow/profiles/production/ "$@"
