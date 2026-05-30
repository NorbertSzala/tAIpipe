rule prepare_metadata_samples:
    input:
        # runs compilation after all summary.tsv for cubar are ready
        summaries = expand(f"{PER_GENOME}/{{sample}}/codon_metrics/{{sample}}_summary.tsv", sample=SAMPLES),
        dataset_sheet = config["paths"]["metadata_dataset"],
        cds_dir = DATA_CDS,
        gcode_dir = GCODES_TRNASCANSE # Pobiera ścieżkę do data/tutorial_data/input/trnascanse
    output:
        # generate file mentioned in DATA_MAINTABLE
        master_table = config["paths"]["metadata_samples"]
    log:
        f"{LOGS}/aggregated/compile_metadata_samples.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/compile_metadata_samples.py \
            --dataset-sheet {input.dataset_sheet} \
            --per-genome-dir {PER_GENOME} \
            --cds-dir {input.cds_dir} \
            --gcode-dir {input.gcode_dir} \
            --output {output.master_table} > {log} 2>&1
        """