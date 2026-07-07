"""
Creates the local Gene Ontology term dictionary used to annotate enrichment results. It uses an existing GO OBO file when available, downloads it only when necessary, and converts the ontology records into a compact tab-separated lookup table.
"""

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
        obo_url = config['go_dictionary'].get("obo_url", "")

    log:
        f"{LOGS}/resources/prepare_go_dictionary.log"

    conda:
        "../envs/python.yaml"

    message:
        "Preparing GO term dictionary"

    shell:
        """
        mkdir -p "$(dirname {log})" "$(dirname {output.tsv}")
        
        python workflow/scripts/resources/parse_go_obo_to_tsv.py \
            --input-obo {params.obo_path} \
            --output-tsv {output.tsv} \
            --obo-url {params.obo_url} \
            > {log} 2>&1
        """