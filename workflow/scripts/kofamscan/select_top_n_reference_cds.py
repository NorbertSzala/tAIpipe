#!/usr/bin/env python3
"""
Select top-N ribosomal CDS sequences per KO from significant KOfam hits.

Input:
    1. ribosome_significant_hits.tsv
       Required columns:
       sample, gene_id, cds_id, protein_id, ko, threshold, score, evalue, definition

    2. CDS FASTA file

Output:
    1. Text file with selected CDS IDs
    2. FASTA file with selected CDS records

Selection logic:
    - group hits by KO
    - sort within each KO by score descending, then evalue ascending
    - keep top N hits per KO
    - deduplicate selected CDS IDs while preserving ranking order
    - extract matching CDS records from input FASTA

This script is intended to be run by Snakemake via the global `snakemake`
object, but it also supports direct command-line execution for debugging.
"""


import argparse
import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class Hit:
    sample: str
    gene_id: str
    cds_id: str
    protein_id: str
    ko: str
    score: float
    evalue: float
    definition: str


@dataclass
class FastaRecord:
    record_id: str
    header: str
    sequence: str


def parse_float(value: str, default: float = float("nan")) -> float:
    """Parse a float value robustly."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_hits(path: Path) -> List[Hit]:
    """Read significant KOfam hits from TSV."""
    required_columns = {
        "sample",
        "gene_id",
        "cds_id",
        "protein_id",
        "ko",
        "score",
        "evalue",
        "definition",
    }

    if not path.exists():
        raise FileNotFoundError(f"Missing hits table: {path}")

    if path.stat().st_size == 0:
        raise ValueError(f"Hits table is empty: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Hits table has no header: {path}")

        missing = required_columns - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Missing required columns in {path}: "
                f"{', '.join(sorted(missing))}"
            )

        hits: List[Hit] = []

        for row in reader:
            cds_id = str(row.get("cds_id", "")).strip()
            ko = str(row.get("ko", "")).strip()

            if not cds_id or not ko:
                continue

            hits.append(
                Hit(
                    sample=str(row.get("sample", "")).strip(),
                    gene_id=str(row.get("gene_id", "")).strip(),
                    cds_id=cds_id,
                    protein_id=str(row.get("protein_id", "")).strip(),
                    ko=ko,
                    score=parse_float(row.get("score", "")),
                    evalue=parse_float(row.get("evalue", "")),
                    definition=str(row.get("definition", "")).strip(),
                )
            )

    if not hits:
        raise ValueError(
            f"No valid hits found in {path}. "
            "The file exists, but no rows with non-empty cds_id and ko were found."
        )

    return hits


def select_top_n_per_ko(hits: List[Hit], top_n: int) -> List[Hit]:
    """Select top-N hits per KO by score desc and evalue asc."""
    if top_n < 1:
        raise ValueError(f"top_n must be >= 1, got: {top_n}")

    grouped: Dict[str, List[Hit]] = defaultdict(list)

    for hit in hits:
        grouped[hit.ko].append(hit)

    selected: List[Hit] = []

    for ko in sorted(grouped):
        ko_hits = grouped[ko]

        ko_hits_sorted = sorted(
            ko_hits,
            key=lambda h: (
                -h.score,
                h.evalue,
                h.cds_id,
            ),
        )

        selected.extend(ko_hits_sorted[:top_n])

    if not selected:
        raise ValueError("No hits selected after top-N filtering.")

    return selected


def deduplicate_hits_by_cds_id(hits: List[Hit]) -> List[Hit]:
    """Deduplicate selected hits by CDS ID while preserving order."""
    seen = set()
    unique_hits: List[Hit] = []

    for hit in hits:
        if hit.cds_id in seen:
            continue

        seen.add(hit.cds_id)
        unique_hits.append(hit)

    return unique_hits


def read_fasta(path: Path) -> List[FastaRecord]:
    """Read FASTA records without external dependencies."""
    if not path.exists():
        raise FileNotFoundError(f"Missing CDS FASTA file: {path}")

    records: List[FastaRecord] = []
    header: Optional[str] = None
    seq_chunks: List[str] = []

    def flush_record() -> None:
        nonlocal header, seq_chunks

        if header is None:
            return

        if not seq_chunks:
            raise ValueError(f"FASTA record has no sequence: {header}")

        record_id = header.split()[0]
        records.append(
            FastaRecord(
                record_id=record_id,
                header=header,
                sequence="".join(seq_chunks),
            )
        )

        header = None
        seq_chunks = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith(">"):
                flush_record()
                header = line[1:].strip()
            else:
                seq_chunks.append(line.strip())

        flush_record()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records


def build_fasta_index(records: List[FastaRecord]) -> Dict[str, FastaRecord]:
    """
    Build a permissive FASTA index.

    The primary key is the first FASTA token. Additional keys are extracted from
    common NCBI-style annotations such as protein_id= and transcript_id=.
    """
    index: Dict[str, FastaRecord] = {}

    for record in records:
        keys = set()

        keys.add(record.record_id)
        keys.add(record.header)

        # Common NCBI-style bracket annotations:
        # [protein_id=KAF9174317.1]
        # [locus_tag=BGX20_000106]
        # [db_xref=...]
        for match in re.finditer(r"\[([A-Za-z0-9_.:-]+)=([^\]]+)\]", record.header):
            keys.add(match.group(2).strip())

        # Also index tokens without punctuation cleanup.
        for token in record.header.split():
            keys.add(token.strip())

        for key in keys:
            if key and key not in index:
                index[key] = record

    return index


def find_record_for_hit(
    hit: Hit,
    fasta_index: Dict[str, FastaRecord],
    records: List[FastaRecord],
) -> Optional[FastaRecord]:
    """
    Match hit to FASTA record.

    Preferred matching:
        1. exact cds_id
        2. exact protein_id
        3. cds_id contained in FASTA header
        4. protein_id contained in FASTA header
    """
    if hit.cds_id in fasta_index:
        return fasta_index[hit.cds_id]

    if hit.protein_id and hit.protein_id in fasta_index:
        return fasta_index[hit.protein_id]

    for record in records:
        if hit.cds_id and hit.cds_id in record.header:
            return record

        if hit.protein_id and hit.protein_id in record.header:
            return record

    return None


def write_ids(path: Path, hits: List[Hit]) -> None:
    """Write selected CDS IDs, one per line."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for hit in hits:
            handle.write(f"{hit.cds_id}\n")


def wrap_sequence(sequence: str, width: int = 80) -> Iterable[str]:
    """Yield wrapped FASTA sequence lines."""
    for i in range(0, len(sequence), width):
        yield sequence[i : i + width]


def write_fasta(path: Path, records: List[Tuple[Hit, FastaRecord]]) -> None:
    """Write selected FASTA records."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for hit, record in records:
            # Keep original FASTA header, but prepend a compact annotation that
            # makes the selected KO/gene traceable.
            header = (
                f"{record.header} "
                f"[selected_ko={hit.ko}] "
                f"[selected_gene_id={hit.gene_id}]"
            )

            handle.write(f">{header}\n")

            for line in wrap_sequence(record.sequence):
                handle.write(f"{line}\n")


def run(
    hits_path: Path,
    cds_fasta_path: Path,
    ids_output_path: Path,
    fasta_output_path: Path,
    top_n: int = 2,
    allow_empty: bool = False,
) -> None:
    """Main execution function."""
    hits = read_hits(hits_path)
    selected_hits = select_top_n_per_ko(hits, top_n=top_n)
    selected_hits = deduplicate_hits_by_cds_id(selected_hits)

    records = read_fasta(cds_fasta_path)
    fasta_index = build_fasta_index(records)

    selected_records: List[Tuple[Hit, FastaRecord]] = []
    missing_hits: List[Hit] = []

    for hit in selected_hits:
        record = find_record_for_hit(hit, fasta_index, records)

        if record is None:
            missing_hits.append(hit)
        else:
            selected_records.append((hit, record))

    if missing_hits:
        preview = "\n".join(
            f"  cds_id={h.cds_id}, protein_id={h.protein_id}, ko={h.ko}"
            for h in missing_hits[:20]
        )

        raise ValueError(
            f"Could not find {len(missing_hits)} selected CDS IDs in FASTA "
            f"{cds_fasta_path}.\n"
            f"First missing hits:\n{preview}"
        )

    if not selected_records:
        if allow_empty:
            ids_output_path.parent.mkdir(parents=True, exist_ok=True)
            fasta_output_path.parent.mkdir(parents=True, exist_ok=True)
            ids_output_path.write_text("", encoding="utf-8")
            fasta_output_path.write_text("", encoding="utf-8")
            print(
                "WARNING: no selected records; empty output files were created.",
                file=sys.stderr,
            )
            return

        raise ValueError(
            "No CDS records selected. Refusing to create empty ribosomal "
            "reference because downstream CAI would be unreliable."
        )

    write_ids(ids_output_path, [hit for hit, _ in selected_records])
    write_fasta(fasta_output_path, selected_records)

    n_ko = len({hit.ko for hit in hits})
    n_selected_ko = len({hit.ko for hit, _ in selected_records})

    print(
        "Selected ribosomal CDS records:",
        f"input_hits={len(hits)}",
        f"input_kos={n_ko}",
        f"top_n={top_n}",
        f"unique_selected_cds={len(selected_records)}",
        f"selected_kos={n_selected_ko}",
        f"ids_output={ids_output_path}",
        f"fasta_output={fasta_output_path}",
        sep="\n  ",
        file=sys.stderr,
    )


def parse_args() -> argparse.Namespace:
    """CLI parser for standalone debugging."""
    parser = argparse.ArgumentParser(
        description="Select top-N ribosomal CDS sequences per KO."
    )

    parser.add_argument(
        "--hits",
        required=True,
        type=Path,
        help="Input ribosome_significant_hits.tsv",
    )

    parser.add_argument(
        "--cds-fasta",
        required=True,
        type=Path,
        help="Input CDS FASTA file",
    )

    parser.add_argument(
        "--ids-output",
        required=True,
        type=Path,
        help="Output selected CDS IDs file",
    )

    parser.add_argument(
        "--fasta-output",
        required=True,
        type=Path,
        help="Output selected CDS FASTA file",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=2,
        help="Number of top hits to keep per KO.",
    )

    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Create empty outputs instead of failing when no records are selected.",
    )

    return parser.parse_args()


def run_from_snakemake() -> None:
    """Run from Snakemake global object."""
    hits_path = Path(snakemake.input["hits"])
    cds_fasta_path = Path(snakemake.input["cds"])

    ids_output_path = Path(snakemake.output["ids"])
    fasta_output_path = Path(snakemake.output["fasta"])

    top_n = int(snakemake.params.get("top_n", 2))
    allow_empty = bool(snakemake.params.get("allow_empty", False))

    run(
        hits_path=hits_path,
        cds_fasta_path=cds_fasta_path,
        ids_output_path=ids_output_path,
        fasta_output_path=fasta_output_path,
        top_n=top_n,
        allow_empty=allow_empty,
    )


def run_from_cli() -> None:
    """Run from command line."""
    args = parse_args()

    run(
        hits_path=args.hits,
        cds_fasta_path=args.cds_fasta,
        ids_output_path=args.ids_output,
        fasta_output_path=args.fasta_output,
        top_n=args.top_n,
        allow_empty=args.allow_empty,
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        run_from_cli()