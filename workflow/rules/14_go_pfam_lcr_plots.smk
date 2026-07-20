# Optional rules for data-driven GO-term plots and PFAM-LCR overlap plots.
# These rules assume that gene_features.tsv already contains GO annotations and
# extended Main_Dataset-derived PFAM/LCR columns.

rule compute_pfam_tai_tail_enrichment:
    input:
        gene_features=config["paths"]["gene_features"],
        pfam_descriptions=config.get("pfam_lcr_plots", {}).get("pfam_description_table", ""),

    output:
        enrichment=f"{DATA_STATISTICS}/pfam_tai_tail_enrichment.tsv",
        qc=f"{DATA_STATISTICS}/pfam_tai_tail_enrichment_qc.tsv",

    params:
        tail_fractions=",".join(map(str, config.get("statistics", {}).get("tai_extremes", {}).get("tail_fractions", [0.10, 0.01]))),
        min_genes_per_sample=config.get("statistics", {}).get("tai_extremes", {}).get("min_genes_per_sample", 100),
        min_tail_size=config.get("statistics", {}).get("tai_extremes", {}).get("min_tail_size", 20),
        min_informative_genomes=config.get("pfam_lcr_plots", {}).get("min_informative_genomes", 3),
        min_total_genes_with_term=config.get("pfam_lcr_plots", {}).get("min_total_genes_with_term", 10),
        fdr_method=config.get("statistics", {}).get("fdr_method", "BH"),

    log:
        f"{LOGS}/statistics/compute_pfam_tai_tail_enrichment.log"

    conda:
        "../envs/r_statistics.yaml"

    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.enrichment:q})" "$(dirname {log:q})"

        Rscript workflow/scripts/statistics/compute_pfam_tai_tail_enrichment.R \
          --gene-features {input.gene_features:q} \
          --pfam-description-table {input.pfam_descriptions:q} \
          --tail-fractions {params.tail_fractions:q} \
          --min-genes-per-sample {params.min_genes_per_sample} \
          --min-tail-size {params.min_tail_size} \
          --min-informative-genomes {params.min_informative_genomes} \
          --min-total-genes-with-term {params.min_total_genes_with_term} \
          --fdr-method {params.fdr_method:q} \
          --output {output.enrichment:q} \
          --qc-output {output.qc:q} \
          > {log:q} 2>&1
        """

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
        manifest=f"{DATA_PLOTS}/go_chosen_terms/plot_manifest.tsv",

    params:
        output_dir=f"{DATA_PLOTS}/go_chosen_terms",
        formats=",".join(config.get("plots", {}).get("output_formats", ["png", "pdf"])),
        max_terms=config.get("go_plots", {}).get("max_plot_terms", 24),
        tail_fraction=config.get("go_enrichment", {}).get("tail_fraction", 0.10),
        max_tail_genes=config.get("go_enrichment", {}).get("max_tail_genes", 0),
        min_genes_with_go_per_sample=config.get("go_enrichment", {}).get("min_genes_with_go_per_sample", 10),

    log:
        f"{LOGS}/plots/plot_chosen_go_terms.log"

    conda:
        "../envs/r_plotting.yaml"

    shell:
        """
        set -euo pipefail
        mkdir -p {params.output_dir:q} "$(dirname {log:q})"

        Rscript workflow/scripts/plots/plot_chosen_go_terms.R \
          --chosen-go-terms {input.chosen:q} \
          --gene-table {input.gene_features:q} \
          --max-terms {params.max_terms} \
          --tail-fraction {params.tail_fraction} \
          --max-tail-genes {params.max_tail_genes} \
          --min-genes-with-go-per-sample {params.min_genes_with_go_per_sample} \
          --formats {params.formats:q} \
          --output-dir {params.output_dir:q} \
          --manifest-output {output.manifest:q} \
          > {log:q} 2>&1
        """


rule plot_pfam_lcr_overlap:
    input:
        gene_features=config["paths"]["gene_features"],
        pfam_descriptions=config.get("pfam_lcr_plots", {}).get("pfam_description_table", ""),
        pfam_tail_enrichment=rules.compute_pfam_tai_tail_enrichment.output.enrichment,

    output:
        manifest=f"{DATA_PLOTS}/pfam_lcr_overlap/plot_manifest.tsv",

    params:
        output_dir=f"{DATA_PLOTS}/pfam_lcr_overlap",
        formats=",".join(config.get("plots", {}).get("output_formats", ["png", "pdf"])),
        top_n_pfam=config.get("pfam_lcr_plots", {}).get("top_n_pfam", 25),
        min_group_n=config.get("pfam_lcr_plots", {}).get("min_group_n", 5),

    log:
        f"{LOGS}/plots/plot_pfam_lcr_overlap.log"

    conda:
        "../envs/r_plotting.yaml"

    shell:
        """
        set -euo pipefail
        mkdir -p {params.output_dir:q} "$(dirname {log:q})"

        Rscript workflow/scripts/plots/plot_pfam_lcr_overlap.R \
          --gene-table {input.gene_features:q} \
          --pfam-tail-enrichment {input.pfam_tail_enrichment:q} \
          --top-n-pfam {params.top_n_pfam} \
          --pfam-description-table {input.pfam_descriptions:q} \
          --min-group-n {params.min_group_n} \
          --formats {params.formats:q} \
          --output-dir {params.output_dir:q} \
          --manifest-output {output.manifest:q} \
          > {log:q} 2>&1
        """

