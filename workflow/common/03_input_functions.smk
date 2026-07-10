from glob import glob
from pathlib import Path

def resolve_single_file(pattern):
    matches = glob(pattern)

    if len(matches) == 0:
        raise FileNotFoundError(f"No file found for pattern: {pattern}")
    if len(matches) > 1:
        raise ValueError(f"Multiple files found for pattern: {pattern}")

    return matches[0]

def get_genome(wildcards):
    pattern = samples_df.loc[wildcards.sample, "genome"]
    return resolve_single_file(f"{DATA_GENOME}/{pattern}")

def get_cds(wildcards):
    pattern = samples_df.loc[wildcards.sample, "cds"]
    return resolve_single_file(f"{DATA_CDS}/{pattern}")

def get_proteome(wildcards):
    pattern = samples_df.loc[wildcards.sample, "proteome"]
    return resolve_single_file(f"{DATA_PROTEOME}/{pattern}")

def get_genetic_code(wildcards):
    return int(samples_df.loc[wildcards.sample, "genetic_code"])

def get_trnascanse_code_arg(wildcards):
    """Return the optional tRNAscan-SE genetic-code argument."""
    genetic_code = get_genetic_code(wildcards)

    if genetic_code == 1:
        return ""

    gcode_file = Path(GCODES_TRNASCANSE) / f"{genetic_code}.gcode"

    if not gcode_file.exists():
        raise FileNotFoundError(
            f"Missing tRNAscan-SE genetic-code file: {gcode_file}. "
            "Use a valid transl_table or add the matching .gcode file."
        )

    return f"-g {gcode_file}"

def get_trnascanse_domain_arg(wildcards):
    """Extract and normalize domain from config/sample.tsv."""
    domain = str(samples_df.loc[wildcards.sample, "domain"]).strip().lower()

    domain_map = {
        "eukarya": "-E",
        "eukaryota": "-E",
        "bacteria": "-B",
        "archaea": "-A",
    }

    if domain not in domain_map:
        raise ValueError(
            f"Invalid domain for sample {wildcards.sample}: {domain}. "
            "Allowed values: Eukarya, Bacteria, Archaea."
        )

    return domain_map[domain]

def get_domain(wildcards):
    """Extract and normalize domain from config/sample.tsv."""
    domain = str(samples_df.loc[wildcards.sample, "domain"]).strip().lower()

    domain_map = {
        "eukarya": "Eukarya",
        "eukaryota": "Eukarya",
        "bacteria": "Bacteria",
        "archaea": "Archaea",
    }

    if domain not in domain_map:
        raise ValueError(
            f"Invalid domain for sample {wildcards.sample}: {domain}. "
            "Allowed values: Eukarya, Bacteria, Archaea."
        )

    return domain_map[domain]

def get_go_dictionary_input(wildcards):
    if config.get("go_dictionary", {}).get("enabled", False):
        return config["go_dictionary"]["path"]
    return []



# --- Targets ---
def build_plot_targets():
    if not config.get("plots", {}).get("enabled", True):
        return []

    plot_groups = config.get("plots", {}).get("groups", {})
    targets = []

    if plot_groups.get("tai_extremes", True):
        targets += [
            f"{DATA_PLOTS}/tai_extremes/binary_features.png",
            f"{DATA_PLOTS}/tai_extremes/continuous_features.png",
        ]

    if plot_groups.get("gene_features", True):
        targets += [
            f"{DATA_PLOTS}/gene_features/tai_by_binary_features.png",
            f"{DATA_PLOTS}/gene_features/tai_vs_continuous_features.png",
            f"{DATA_PLOTS}/gene_features/gene_feature_effects.png",
        ]

    if plot_groups.get("genome_metrics", True):
        targets += [
            f"{DATA_PLOTS}/genome_metrics/genome_metrics_by_phylum.png",
            f"{DATA_PLOTS}/genome_metrics/genome_metrics_by_lifestyle.png",
        ]

    if plot_groups.get("codon_profiles", True):
        targets += [
            f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap.png",
        ]

    if plot_groups.get("correlations", False):
        targets += [
            f"{DATA_PLOTS}/correlations/gene_level_correlation_matrix.png",
            f"{DATA_PLOTS}/correlations/genome_level_correlation_matrix.png",
        ]

    if plot_groups.get("qc", False):
        targets += [
            f"{DATA_PLOTS}/qc/trna_qc_flags.png",
        ]

    if (
        plot_groups.get("go_enrichment", False)
        and config.get("go_enrichment", {}).get("enabled", False)
    ):
        targets += [
            f"{DATA_PLOTS}/go_enrichment/go_top_terms_dotplot.png",
            f"{DATA_PLOTS}/go_enrichment/go_odds_ratio_forest.png",
        ]

    return targets

def build_other_targets():
    targets = []

    if (
        config.get("go_enrichment", {}).get("enabled", False)
        and config.get("go_plots", {}).get("enabled", False)
    ):
        targets += [
            f"{DATA_PLOTS}/go_chosen_terms",
        ]

    if config.get("pfam_lcr_plots", {}).get("enabled", False):
        targets += [
            f"{DATA_PLOTS}/pfam_lcr_overlap",
        ]

    return targets