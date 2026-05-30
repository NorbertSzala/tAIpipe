rule aggregate_and_report:
    input:
        # Ensures that the report is triggered only when all summary.tsv files are ready
        summaries = expand(f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv", sample=SAMPLES),
        # Pass specific paths needed for the R script
        per_genome_dir = PER_GENOME,
        # Safely read the corrected path directly from the configuration object
        metadata_master = config["paths"]["metadata_samples"],
        samples_sheet = config["paths"]["metadata_dataset"],
        # IMPORTANT: to ensure updates when report structure is changes
        template = "workflow/scripts/report_template.Rmd"
    output:
        html_report = f"{REPORT_FILE}/summary_report.html",
        md_report = f"{REPORT_FILE}/summary_report.md"
    log:
        f"{LOGS}/aggregated/generate_report.log"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/generate_report.R"