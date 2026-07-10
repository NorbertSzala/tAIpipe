# Optional rules for data-driven GO-term plots and PFAM-LCR overlap plots.
# These rules assume that gene_features.tsv already contains GO annotations and
# extended Main_Dataset-derived PFAM/LCR columns.

rule select_chosen_go_terms:
    input:
        go_enrichment=config["paths"]["go_enrichment"],
        gene_features=config["paths"]["gene_features"],
        go_dictionary=config["go_dictionary"]["path"],

    output:
        chosen=f"{DATA_STATISTICS}/chosen_GOterms.tsv",
        diagnostics=f"{DATA_STATISTICS}/chosen_GOterms_diagnostics.tsv",

    params:
        q_threshold=config.get("go_plots", {}).get("q_threshold", 0.05),
        min_abs_log2_or=config.get("go_plots", {}).get("min_abs_log2_or", 0.25),
        min_informative_genomes=config.get("go_plots", {}).get("min_informative_genomes", 5),
        min_total_genes_with_term=config.get("go_plots", {}).get("min_total_genes_with_term", 30),
        max_terms_per_tail_namespace=config.get("go_plots", {}).get("max_terms_per_tail_namespace", 8),
        namespaces=",".join(config.get("go_plots", {}).get("namespaces", [
            "biological_process",
            "molecular_function",
            "cellular_component",
        ])),

    log:
        f"{LOGS}/statistics/select_chosen_go_terms.log"

    conda:
        "../envs/r_statistics.yaml"

    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.chosen:q})" "$(dirname {log:q})"

        Rscript workflow/scripts/statistics/select_chosen_go_terms.R \
          --go-enrichment {input.go_enrichment:q} \
          --gene-features {input.gene_features:q} \
          --go-dictionary {input.go_dictionary:q} \
          --q-threshold {params.q_threshold} \
          --min-abs-log2-or {params.min_abs_log2_or} \
          --min-informative-genomes {params.min_informative_genomes} \
          --min-total-genes-with-term {params.min_total_genes_with_term} \
          --max-terms-per-tail-namespace {params.max_terms_per_tail_namespace} \
          --namespaces {params.namespaces:q} \
          --output {output.chosen:q} \
          --diagnostics-output {output.diagnostics:q} \
          > {log:q} 2>&1
        """


rule plot_chosen_go_terms:
    input:
        chosen=rules.select_chosen_go_terms.output.chosen,
        gene_features=config["paths"]["gene_features"],

    output:
        directory(f"{DATA_PLOTS}/go_chosen_terms"),

    params:
        max_terms=config.get("go_plots", {}).get("max_plot_terms", 24),

    log:
        f"{LOGS}/plots/plot_chosen_go_terms.log"

    conda:
        "../envs/r_plotting.yaml"

    shell:
        """
        set -euo pipefail
        mkdir -p {output:q} "$(dirname {log:q})"

        Rscript workflow/scripts/plots/plot_chosen_go_terms.R \
          --chosen-go-terms {input.chosen:q} \
          --gene-table {input.gene_features:q} \
          --max-terms {params.max_terms} \
          --output-dir {output:q} \
          > {log:q} 2>&1
        """


rule plot_pfam_lcr_overlap:
    input:
        gene_features=config["paths"]["gene_features"],

    output:
        directory(f"{DATA_PLOTS}/pfam_lcr_overlap"),

    params:
        top_n_pfam=config.get("pfam_lcr_plots", {}).get("top_n_pfam", 25),

    log:
        f"{LOGS}/plots/plot_pfam_lcr_overlap.log"

    conda:
        "../envs/r_plotting.yaml"

    shell:
        """
        set -euo pipefail
        mkdir -p {output:q} "$(dirname {log:q})"

        Rscript workflow/scripts/plots/plot_pfam_lcr_overlap.R \
          --gene-table {input.gene_features:q} \
          --top-n-pfam {params.top_n_pfam} \
          --output-dir {output:q} \
          > {log:q} 2>&1
        """
