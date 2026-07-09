from glob import glob
from pathlib import Path

def resolve_single_file(pattern):
    matches = glob(pattern)

    if len(matches) == 0:
        raise FileNotFoundError(f"No file found for pattern: {pattern}")
    if len(matches) > 1:
        raise ValueError(f"Multiple files found for pattern: {pattern}")

    return matches[0]

def get_genome(wildcards):
    pattern = samples_df.loc[wildcards.sample, "genome"]
    return resolve_single_file(f"{DATA_GENOME}/{pattern}")

def get_cds(wildcards):
    pattern = samples_df.loc[wildcards.sample, "cds"]
    return resolve_single_file(f"{DATA_CDS}/{pattern}")

def get_proteome(wildcards):
    pattern = samples_df.loc[wildcards.sample, "proteome"]
    return resolve_single_file(f"{DATA_PROTEOME}/{pattern}")

def get_genetic_code(wildcards):
    return int(samples_df.loc[wildcards.sample, "genetic_code"])

def get_trnascanse_code_arg(wildcards):
    """Return the optional tRNAscan-SE genetic-code argument."""
    genetic_code = get_genetic_code(wildcards)

    if genetic_code == 1:
        return ""

    gcode_file = Path(GCODES_TRNASCANSE) / f"{genetic_code}.gcode"

    if not gcode_file.exists():
        raise FileNotFoundError(
            f"Missing tRNAscan-SE genetic-code file: {gcode_file}. "
            "Use a valid transl_table or add the matching .gcode file."
        )

    return f"-g {gcode_file}"

def get_domain(wildcards):
    """Extract and normalize domain from metadata_dataset.tsv."""
    domain = str(samples_df.loc[wildcards.sample, "domain"]).strip().lower()

    domain_map = {
        "eukarya": "Eukarya",
        "eukaryota": "Eukarya",
        "bacteria": "Bacteria",
        "archaea": "Archaea",
    }

    if domain not in domain_map:
        raise ValueError(
            f"Invalid domain for sample {wildcards.sample}: {domain}. "
            "Allowed values: Eukarya, Bacteria, Archaea."
        )

    return domain_map[domain]

def get_go_dictionary_input(wildcards):
    if config.get("go_dictionary", {}).get("enabled", False):
        return config["go_dictionary"]["path"]
    return []
