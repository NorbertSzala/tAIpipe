rule prepare_metadata_samples:
    input:
        # runs compilation after all summary.tsv for cubar are ready
        summaries = expand(f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv", sample=SAMPLES),
        metadata_dataset = config["paths"]["metadata_dataset"],
        cds_dir = DATA_CDS,
        gcode_dir = GCODES_TRNASCANSE # Pobiera ścieżkę do data/tutorial_data/input/trnascanse
    output:
        # generate file mentioned in DATA_MAINTABLE
        metadata_samples = config["paths"]["metadata_samples"]
    log:
        f"{LOGS}/aggregated/compile_metadata_samples.log"
    conda:
        "../envs/r.yaml"
    shell:
        """
        Rscript workflow/scripts/compile_metadata_samples.R \
            --metadata_dataset {input.metadata_dataset} \
            --per-genome-dir {PER_GENOME} \
            --cds-dir {input.cds_dir} \
            --gcode-dir {input.gcode_dir} \
            --output {output.metadata_samples} > {log} 2>&1
        """