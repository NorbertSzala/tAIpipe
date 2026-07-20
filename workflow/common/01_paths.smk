from pathlib import Path

# Path aliases used across rules. Keep this file small and declarative.
DATA_CDS = config["paths"]["data_cds"]
DATA_GENOME = config["paths"]["data_genome"]
DATA_PROTEOME = config["paths"]["data_proteome"]
GCODES_TRNASCANSE = config["paths"]["genetic_codes_trnascanse"]
GCODES_TABLES = config["paths"].get("genetic_codes_tables")

PER_GENOME = config["paths"]["per_genome"]
LOGS = config["paths"]["logs"]
BENCHMARKS = config["paths"]["benchmarks"]

DATA_PLOTS = config["paths"]["plots"]
DATA_TABLES = config["paths"]["tables"]
DATA_QC = config["paths"]["qc"]
DATA_STATISTICS = config["paths"]["statistics"]
REPORT_FILE = config["paths"]["reports"]

GENE_FEATURES = config["paths"]["gene_features"]
GENOME_SUMMARY = config["paths"]["genome_summary"]
CODON_PROFILES = config["paths"]["codon_profiles"]

TEMPLATE_PATH = config.get("report_config", {}).get("template_path")

def require_existing_path(path, label):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Required {label} does not exist: {path}")
    return path

require_existing_path(DATA_CDS, "CDS directory")
require_existing_path(DATA_GENOME, "genome directory")
require_existing_path(DATA_PROTEOME, "proteome directory")
require_existing_path(GCODES_TRNASCANSE, "tRNAscan-SE genetic codes directory")
require_existing_path(config["paths"]["metadata_dataset"], "metadata dataset file")
