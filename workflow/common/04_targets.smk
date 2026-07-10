# Target lists only. Do not use directory() in rule all inputs.
# Keep these paths synchronized with rule outputs in workflow/rules/*.smk.

GO_ENRICHMENT_OUTPUTS = []
if config.get("go_enrichment", {}).get("enabled", False):
    GO_ENRICHMENT_OUTPUTS = [
        config["paths"]["go_enrichment"],
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

REFERENCE_SELECTION = config.get("ribosomal_reference", {}).get(
    "selection_strategy", "best_per_ko"
)
REFERENCE_TOP_N = int(config.get("ribosomal_reference", {}).get("top_n_per_ko", 1))

if REFERENCE_SELECTION == "best_per_ko" or REFERENCE_TOP_N == 1:
    REFERENCE_SUFFIX = "best_per_ko"
elif REFERENCE_SELECTION == "top_n_per_ko":
    REFERENCE_SUFFIX = f"top{REFERENCE_TOP_N}_per_ko"
else:
    raise ValueError(
        "ribosomal_reference.selection_strategy must be best_per_ko or top_n_per_ko"
    )

KOFAM_TARGETS = expand(
    f"{PER_GENOME}/{{sample}}/kofamscan/ribosomal_reference_cds.{REFERENCE_SUFFIX}.fna",
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
    + expand(
        f"{PER_GENOME}/{{sample}}/qc/{{sample}}_metric_qc.tsv",
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

if config.get("statistics", {}).get("lifestyle_within_phylum_permutations", {}).get("enabled", False):
    STATISTICS_TARGETS += [
        f"{DATA_STATISTICS}/lifestyle_within_phylum_permutation_tests.tsv",
    ]

if config.get("statistics", {}).get("per_genome_feature_effects", {}).get("enabled", False):
    STATISTICS_TARGETS += [
        f"{DATA_STATISTICS}/gene_feature_per_genome_effects.tsv",
        f"{DATA_STATISTICS}/gene_feature_effect_meta_tests.tsv",
    ]


def build_plot_targets():
    if not config.get("plots", {}).get("enabled", True):
        return []

    plot_groups = config.get("plots", {}).get("groups", {})
    formats = config.get("plots", {}).get("output_formats", ["png", "pdf"])
    targets = []

    def add_plot(stem):
        for fmt in formats:
            targets.append(f"{stem}.{fmt}")

    if plot_groups.get("tai_extremes", True):
        add_plot(f"{DATA_PLOTS}/tai_extremes/binary_features")
        add_plot(f"{DATA_PLOTS}/tai_extremes/continuous_features")

    if plot_groups.get("gene_features", True):
        add_plot(f"{DATA_PLOTS}/gene_features/tai_by_binary_features")
        add_plot(f"{DATA_PLOTS}/gene_features/tai_vs_continuous_features")
        add_plot(f"{DATA_PLOTS}/gene_features/gene_feature_effects")

    if plot_groups.get("genome_metrics", True):
        add_plot(f"{DATA_PLOTS}/genome_metrics/genome_metrics_by_phylum")
        add_plot(f"{DATA_PLOTS}/genome_metrics/genome_metrics_by_lifestyle")

    if plot_groups.get("codon_profiles", True):
        min_group_n = int(config.get("codon_profile_plots", {}).get("heatmap_min_group_n", 5))
        add_plot(f"{DATA_PLOTS}/codon_profiles/codon_usage_variability")
        add_plot(f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_phylum_large_codon_x")
        add_plot(f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_lifestyle_large_codon_x")
        add_plot(f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_phylum_count{min_group_n}_large_codon_x")
        add_plot(f"{DATA_PLOTS}/codon_profiles/trna_weights_heatmap_by_lifestyle_count{min_group_n}_large_codon_x")
        targets.append(f"{DATA_PLOTS}/codon_profiles/codon_usage_variability_method.tsv")

    if plot_groups.get("correlations", False):
        add_plot(f"{DATA_PLOTS}/correlations/gene_level_correlation_matrix")
        add_plot(f"{DATA_PLOTS}/correlations/genome_level_correlation_matrix")

    if plot_groups.get("qc", False):
        add_plot(f"{DATA_PLOTS}/qc/trna_qc_flags")

    if (
        plot_groups.get("go_enrichment", False)
        and config.get("go_enrichment", {}).get("enabled", False)
    ):
        add_plot(f"{DATA_PLOTS}/go_enrichment/go_top_terms_dotplot")
        add_plot(f"{DATA_PLOTS}/go_enrichment/go_odds_ratio_forest")

    return targets


def build_script_suggestion_plot_targets():
    if not config.get("script_suggestion_plots", {}).get("enabled", False):
        return []

    targets = [
        f"{DATA_PLOTS}/script_suggestions/features_tai",
        f"{DATA_PLOTS}/script_suggestions/go_term_features",
        f"{DATA_PLOTS}/script_suggestions/codon_trna",
        f"{DATA_PLOTS}/script_suggestions/genome_gc_tai",
        f"{DATA_PLOTS}/script_suggestions/correlations",
    ]

    return targets


def build_other_targets():
    targets = []

    if (
        config.get("go_enrichment", {}).get("enabled", False)
        and config.get("go_plots", {}).get("enabled", False)
    ):
        targets += [f"{DATA_PLOTS}/go_chosen_terms"]

    if config.get("pfam_lcr_plots", {}).get("enabled", False):
        targets += [f"{DATA_PLOTS}/pfam_lcr_overlap"]

    return targets


PLOT_TARGETS = build_plot_targets()
SCRIPT_SUGGESTION_PLOT_TARGETS = build_script_suggestion_plot_targets()
OTHER_TARGETS = build_other_targets()

WORKFLOW_TARGETS = (
    PER_SAMPLE_TARGETS
    + KOFAM_TARGETS
    + CANONICAL_TABLE_TARGETS
    + STATISTICS_TARGETS
    + GO_ENRICHMENT_OUTPUTS
    + PLOT_TARGETS
    + SCRIPT_SUGGESTION_PLOT_TARGETS
    + OTHER_TARGETS
)
