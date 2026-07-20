#!/usr/bin/env bash
set -Eeuo pipefail

DRY_RUN=0
RUN_PHYLOGENY=0
CLEAN_OUTPUTS=1
CORES="${CORES:-32}"

usage() {
  cat <<'EOF'
Usage:
  bash rerun_changed_outputs.sh [--dry-run] [--with-phylogeny] [--no-clean]

Options:
  --dry-run          Show the Snakemake execution plan without running jobs.
  --with-phylogeny   Also rebuild KO prevalence and the ribosomal species tree.
  --no-clean         Do not remove previous plot directories before rerunning.
                     Declared outputs are still forced by Snakemake.

Environment:
  CORES=32           Number of cores used by Snakemake and phylogeny tools.
EOF
}

while (($#)); do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      ;;
    --with-phylogeny)
      RUN_PHYLOGENY=1
      ;;
    --no-clean)
      CLEAN_OUTPUTS=0
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

if [[ ! -f workflow/Snakefile ]]; then
  echo "Run this script from the tAIpipe repository root." >&2
  exit 1
fi

PLOT_RULES=(
  plot_tai_extreme_features
  plot_gene_feature_overview
  plot_genome_metric_overview
  plot_codon_profile_overview
  plot_go_enrichment_overview
  plot_correlation_overview
  plot_qc_overview
  plot_script_suggestion_feature_tai
  plot_script_suggestion_codon_trna
  plot_script_suggestion_genome_gc_tai
  plot_script_suggestion_correlations
  plot_script_suggestion_go_term_features
  select_chosen_go_terms
  plot_chosen_go_terms
  plot_pfam_lcr_overlap
)

CHANGED_OUTPUT_DIRS=(
  results/plots/tai_extremes
  results/plots/gene_features
  results/plots/genome_metrics
  results/plots/codon_profiles
  results/plots/go_enrichment
  results/plots/correlations
  results/plots/qc
  results/plots/go_chosen_terms
  results/plots/pfam_lcr_overlap
  results/plots/script_suggestions/features_tai
  results/plots/script_suggestions/codon_trna
  results/plots/script_suggestions/genome_gc_tai
  results/plots/script_suggestions/correlations
  results/plots/script_suggestions/go_term_features
)

if ((DRY_RUN == 0 && CLEAN_OUTPUTS == 1)); then
  echo "Removing only outputs produced by changed plotting scripts..."
  rm -rf -- "${CHANGED_OUTPUT_DIRS[@]}"
fi

SNAKEMAKE_ARGS=(
  -s workflow/Snakefile
  --profile workflow/profiles/production
  --cores "$CORES"
  --rerun-incomplete
  --printshellcmds
  "${PLOT_RULES[@]}"
  --forcerun "${PLOT_RULES[@]}"
)

if ((DRY_RUN == 1)); then
  SNAKEMAKE_ARGS+=(--dry-run)
fi

echo "Running changed plotting rules..."
snakemake "${SNAKEMAKE_ARGS[@]}"

if ((DRY_RUN == 1)); then
  exit 0
fi

echo "Rebuilding ribosomal KO prevalence summary..."
python workflow/scripts/resources/summarize_ribosomal_ko_prevalence.py \
  --samples config/samples.tsv \
  --results-root results/per_genome \
  --output results/tables/ribosomal_ko_prevalence.tsv

if ((RUN_PHYLOGENY == 0)); then
  echo
  echo "Plots and KO prevalence rebuilt."
  echo "Use --with-phylogeny to also rebuild results/phylogeny/ribosomal_markers."
  exit 0
fi

ENV_NAME="taipipe-ribosomal-phylogeny"
ENV_FILE="workflow/envs/ribosomal_phylogeny.yaml"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda is required for --with-phylogeny." >&2
  exit 1
fi

if conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "$ENV_NAME"; then
  echo "Updating Conda environment: $ENV_NAME"
  if command -v mamba >/dev/null 2>&1; then
    mamba env update -n "$ENV_NAME" -f "$ENV_FILE" --prune
  else
    conda env update -n "$ENV_NAME" -f "$ENV_FILE" --prune
  fi
else
  echo "Creating Conda environment: $ENV_NAME"
  if command -v mamba >/dev/null 2>&1; then
    mamba env create -f "$ENV_FILE"
  else
    conda env create -f "$ENV_FILE"
  fi
fi

echo "Removing previous phylogeny output..."
rm -rf -- results/phylogeny/ribosomal_markers

echo "Building ribosomal marker species tree..."
conda run -n "$ENV_NAME" \
  python workflow/scripts/phylogeny/build_ribosomal_phylogeny.py \
    --samples config/samples.tsv \
    --results-root results/per_genome \
    --output-dir results/phylogeny/ribosomal_markers \
    --marker-prevalence 0.90 \
    --final-marker-coverage 0.80 \
    --min-genome-marker-coverage 0.70 \
    --threads "$CORES" \
    --ufboot 1000 \
    --alrt 1000

echo
echo "All changed outputs rebuilt."
echo "Tree: results/phylogeny/ribosomal_markers/fungal_ribosomal_species_tree.treefile"