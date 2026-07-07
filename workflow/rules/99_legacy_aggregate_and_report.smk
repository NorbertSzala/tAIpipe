"""
Legacy reporting rule that collects per-genome outputs and renders the combined HTML and Markdown reports. It depends on the older aggregation layout and should remain disabled until it is updated to consume the canonical tables and the current report template.
"""

rule aggregate_and_report:
    input:
        summaries=expand(
            f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv",
            sample=SAMPLES,
        ),
        rscu=expand(
            f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_rscu.csv",
            sample=SAMPLES,
        ),
        trna_weights=expand(
            f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_trna_weights.csv",
            sample=SAMPLES,
        ),
        samples=config["paths"]["metadata_samples"],
        dataset=config["paths"]["metadata_dataset"],
        template=TEMPLATE_PATH,
    output:
        html_report=f"{REPORT_FILE}/summary_report.html",
        md_report=f"{REPORT_FILE}/summary_report.md",
        assets=directory(f"{REPORT_FILE}/summary_report_files"),
    log:
        f"{LOGS}/aggregated/generate_report.log"
    params:
        per_genome_dir=PER_GENOME,
        report_file_path=REPORT_FILE,
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/report/generate_report.R"
