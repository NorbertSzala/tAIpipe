def get_obo_url_arg(wildcards):
    obo_url = config.get("go_dictionary", {}).get("obo_url", "")
    if obo_url:
        return f'--obo-url "{obo_url}"'
    return ""


rule prepare_go_dictionary:
    """
    Prepare GO term dictionary.

    If TSV exists, the script validates it and exits.
    If TSV is missing, the script checks for local OBO.
    If OBO is missing, it downloads OBO from the configured URL.
    Then it converts OBO to TSV.
    """
    output:
        tsv=config['go_dictionary']['path']

    params:
        obo_path = config['go_dictionary']['obo_path'],
        obo_url_arg = get_obo_url_arg

    log:
        f"{LOGS}/resources/prepare_go_dictionary.log"

    conda:
        "../envs/python.yaml"

    message:
        "Preparing GO term dictionary"

    shell:
        """
        mkdir -p $(dirname {log})
        
        python workflow/scripts/parse_go_obo_to_tsv.py \
            --input-obo {params.obo_path} \
            --output-tsv {output.tsv} \
            {params.obo_url_arg} \
            > {log} 2>&1
        """
