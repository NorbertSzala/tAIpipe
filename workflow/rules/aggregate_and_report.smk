rule aggregate_and_report:
    """
    Create interactive raport saved in HTML and .Rmd
    """
    
    input:
        # Ensures that the report is triggered only when all summary.tsv files are ready
        summaries = expand(f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv", sample=SAMPLES),
        # Pass specific paths needed for the R script
        # Safely read the corrected path directly from the configuration object
        samples = config["paths"]["metadata_samples"],
        dataset = config["paths"]["metadata_dataset"],
        # IMPORTANT: to ensure updates when report structure is changes
        template=TEMPLATE_PATH,
        go_dictionary  = get_go_dictionary_input
    
    output:
        html_report = f"{REPORT_FILE}/summary_report.html",
        md_report = f"{REPORT_FILE}/summary_report.md",
        # assets = directory(f"{REPORT_FILE}/summary_report_files")
    
    params:
        per_genome_dir = PER_GENOME,

    log:
        f"{LOGS}/aggregated/generate_report.log"
    
    params:
        # Transfer static environment configurations and directory locations safely to params
        per_genome_dir = PER_GENOME,
        report_file_path = REPORT_FILE,
        go_dictionary_enabled=lambda wildcards: config.get("go_dictionary", {}).get("enabled", False),
        go_dictionary_path=lambda wildcards: config.get("go_dictionary", {}).get("path", "")
    
    conda:
        "../envs/r.yaml"
    
    script:
        "../scripts/generate_report.R"

