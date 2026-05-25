#!/usr/bin/env bash
set -euo pipefail

Rscript -e 'if (!requireNamespace("cubar", quietly = TRUE)) install.packages("cubar", repos = "https://cloud.r-project.org")'