rule aggregate_and_report:
    input:
        # Ensures that the report is triggered only when all summary.tsv files are ready
        summaries = expand(f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv", sample=SAMPLES),
        # Pass specific paths needed for the R script
        per_genome_dir = PER_GENOME,
        # Safely read the corrected path directly from the configuration object
        metadata_master = config["paths"]["main_dataset"],
        samples_sheet = config["paths"]["samples"]
    output:
        html_report = f"{REPORT_FILE}/summary_report.html",
        md_report = f"{REPORT_FILE}/summary_report.md"
    log:
        f"{LOGS}/aggregated/generate_report.log"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/generate_report.R"