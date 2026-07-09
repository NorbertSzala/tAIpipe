# Optional output groups. This file should contain only target lists.

GO_ENRICHMENT_OUTPUTS = []
if config.get("go_enrichment", {}).get("enabled", False):
    GO_ENRICHMENT_OUTPUTS = [
        config["paths"]["go_enrichment"]
    ]

KOFAM_ENABLED = config.get("kofamscan", {}).get("enabled", False)
if not KOFAM_ENABLED:
    raise ValueError(
        "kofamscan.enabled must be true because codon_usage_metrics requires "
        "a ribosomal reference CDS set for CAI calculation."
    )

require_existing_path(config["kofamscan"]["database"]["ko_list"], "KOfam ko_list")
require_existing_path(config["kofamscan"]["database"]["profiles_dir"], "KOfam profiles directory")
require_existing_path(
    config["kofamscan"]["ribosome_reference"]["brite_json"],
    "KEGG BRITE ko03011 JSON",
)

KOFAM_TARGETS = expand(
    f"{PER_GENOME}/{{sample}}/kofamscan/ribosomal_reference_cds.fna",
    sample=SAMPLES,
)

PER_SAMPLE_TARGETS = (
    expand(
        f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv",
        sample=SAMPLES,
    )
    + expand(
        f"{PER_GENOME}/{{sample}}/qc/{{sample}}_trna_profile_qc.tsv",
        sample=SAMPLES,
    )
)

CANONICAL_TABLE_TARGETS = [
    config["paths"]["gene_features"],
    config["paths"]["genome_summary"],
    config["paths"]["codon_profiles"],
]

STATISTICS_TARGETS = [
    config["paths"]["gene_feature_tests"],
    config["paths"]["genome_group_tests"],
]

PLOT_TARGETS = []
if config.get("plots", {}).get("enabled", True):
    plot_groups = config.get("plots", {}).get("groups", {})

    if plot_groups.get("tai_extremes", True):
        PLOT_TARGETS += [
            "results/plots/tai_extremes/binary_features.png",
            "results/plots/tai_extremes/continuous_features.png",
        ]

    if plot_groups.get("gene_features", True):
        PLOT_TARGETS += [
            "results/plots/gene_features/tai_by_binary_features.png",
            "results/plots/gene_features/tai_vs_continuous_features.png",
            "results/plots/gene_features/gene_feature_effects.png",
        ]

    if plot_groups.get("genome_metrics", True):
        PLOT_TARGETS += [
            "results/plots/genome_metrics/genome_metrics_by_phylum.png",
            "results/plots/genome_metrics/genome_metrics_by_lifestyle.png",
        ]

    if plot_groups.get("codon_profiles", True):
        PLOT_TARGETS += [
            "results/plots/codon_profiles/trna_weights_heatmap.png",
            "results/plots/codon_profiles/rscu_heatmap.png",
        ]

    if plot_groups.get("go_enrichment", False) and config.get("go_enrichment", {}).get("enabled", False):
        PLOT_TARGETS += [
            "results/plots/go_enrichment/go_top_terms_dotplot.png",
            "results/plots/go_enrichment/go_odds_ratio_forest.png",
        ]

    if plot_groups.get("correlations", False):
        PLOT_TARGETS += [
            "results/plots/correlations/gene_level_correlation_matrix.png",
            "results/plots/correlations/genome_level_correlation_matrix.png",
        ]

    if plot_groups.get("qc", False):
        PLOT_TARGETS += [
            "results/plots/qc/trna_qc_flags.png",
        ]

WORKFLOW_TARGETS = (
    PER_SAMPLE_TARGETS
    + KOFAM_TARGETS
    + CANONICAL_TABLE_TARGETS
    + STATISTICS_TARGETS
    + GO_ENRICHMENT_OUTPUTS
    + PLOT_TARGETS
)
