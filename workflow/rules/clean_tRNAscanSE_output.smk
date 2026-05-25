rule convert_trnascanse_output_to_tsv:
    """
    Script converting unfriendly tRNAscan-SE output to more functional .tsv file

    ./scripts/convert_trnascanse_output_to_tsv.py -I results/per_genome/Spombe/trnascan/Spombe_trnascan.out -O ./try.tsv --remove-nnn
    """

    input:
        trnascan_out = f"{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.out",

    output:
        clean_trnascan = f"{PER_GENOME}/{{sample}}/trnascan/{{sample}}_trnascan.tsv" 

    log:
        f'{LOGS}/{{sample}}/convert_trnascanse_output_to_tsv.log'

    params:
        remove_nnn = "--remove-nnn" if config.get('trnascanse_clean', {}).get("remove_nnn", True) else "",

        keep_pseudo = "--keep-pseudo" if config.get('trnascanse_clean', {}).get("keep_pseudo", True) else ""

    conda:
        '../envs/python.yaml'

    message:
        "Cleaning tRNAscan-SE output and saving into easy to interpret .tsv file"

    shell:
        """
        mkdir -p $(dirname {output.clean_trnascan}) $(dirname {log})

        python3 workflow/scripts/convert_trnascanse_output_to_tsv.py \
            -I {input.trnascan_out} \
            -O {output.clean_trnascan} \
            {params.remove_nnn} \
            {params.keep_pseudo}
            > {log} 2>&1
        """
