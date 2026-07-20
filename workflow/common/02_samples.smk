import pandas as pd
from snakemake.utils import validate

samples_df = pd.read_csv(config["paths"]["metadata_dataset"], sep="\t")
validate(samples_df, "../schemas/metadata_dataset.schema.yaml")

def normalize_include(value):
    """Normalize supported boolean encodings from the validated sample sheet."""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)

    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "1", "yes", "y"}:
        return True
    if normalized in {"false", "f", "0", "no", "n"}:
        return False
    raise ValueError(f"Unsupported include value: {value!r}")

samples_df["include"] = samples_df["include"].map(normalize_include)

if samples_df["sample"].duplicated().any():
    duplicated = samples_df.loc[samples_df["sample"].duplicated(), "sample"].tolist()
    raise ValueError(f"Duplicated sample identifiers: {duplicated}")

samples_df = samples_df.set_index("sample")
samples_df = samples_df[samples_df["include"]]

SAMPLES = samples_df.index.tolist()

if not SAMPLES:
    raise ValueError("No samples with include=true were found in the metadata table")
