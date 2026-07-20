#!/usr/bin/env python3
"""Build a low-cost fungal species tree from genome-specific Kofam hits.

The script selects prevalent ribosomal KO markers, keeps one unambiguous best
protein hit per KO and genome, filters obvious length outliers, aligns each
marker with MAFFT, trims alignments with trimAl, concatenates the retained
markers and estimates a partitioned maximum-likelihood tree with IQ-TREE 2.

It consumes outputs already produced by the KofamScan module and does not rerun
homology searches.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

TRUE_VALUES = {"1", "true", "yes", "y", "include"}


@dataclass(frozen=True)
class Sample:
    sample: str
    proteome_pattern: str


@dataclass(frozen=True)
class Hit:
    sample: str
    ko: str
    protein_id: str
    score: float
    threshold: float
    evalue: float
    definition: str

    @property
    def normalized_score(self) -> float:
        if math.isfinite(self.threshold) and self.threshold > 0:
            return self.score / self.threshold
        return self.score


@dataclass(frozen=True)
class SelectedSequence:
    sample: str
    ko: str
    protein_id: str
    definition: str
    score: float
    threshold: float
    normalized_score: float
    second_normalized_score: float | None
    ambiguity_ratio: float | None
    sequence: str


def parse_float(value: object, default: float = math.nan) -> float:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def is_included(value: str | None) -> bool:
    if value is None or not value.strip():
        return True
    return value.strip().lower() in TRUE_VALUES


def read_samples(path: Path) -> list[Sample]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "sample" not in reader.fieldnames:
            raise ValueError(f"{path} must contain a sample column")
        result = []
        for row in reader:
            if not is_included(row.get("include")):
                continue
            sample = str(row.get("sample", "")).strip()
            if sample:
                result.append(Sample(sample=sample, proteome_pattern=str(row.get("proteome", "")).strip()))
    if not result:
        raise ValueError(f"No included samples found in {path}")
    duplicates = [name for name, n in Counter(x.sample for x in result).items() if n > 1]
    if duplicates:
        raise ValueError("Duplicated sample IDs: " + ", ".join(sorted(duplicates)))
    return result


def read_hits(path: Path, expected_sample: str) -> list[Hit]:
    if not path.exists():
        raise FileNotFoundError(f"Missing Kofam hit table: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ko", "protein_id", "score", "evalue"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"{path} lacks columns: {', '.join(missing)}")
        result = []
        for row in reader:
            ko = str(row.get("ko", "")).strip()
            protein_id = str(row.get("protein_id", "")).strip()
            if ko and protein_id:
                result.append(
                    Hit(
                        sample=str(row.get("sample", expected_sample)).strip() or expected_sample,
                        ko=ko,
                        protein_id=protein_id,
                        score=parse_float(row.get("score"), -math.inf),
                        threshold=parse_float(row.get("threshold")),
                        evalue=parse_float(row.get("evalue"), math.inf),
                        definition=str(row.get("definition", "")).strip(),
                    )
                )
    return result


def resolve_proteome(sample: Sample, protein_template: str, proteome_root: Path | None) -> Path:
    rendered = Path(protein_template.format(sample=sample.sample))
    if rendered.exists():
        return rendered
    pattern = sample.proteome_pattern
    if pattern:
        candidate_pattern = Path(pattern)
        if not candidate_pattern.is_absolute() and proteome_root is not None:
            candidate_pattern = proteome_root / candidate_pattern
        matches = sorted(candidate_pattern.parent.glob(candidate_pattern.name))
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError(
                f"Proteome glob for {sample.sample} matched multiple files: "
                + ", ".join(str(x) for x in matches[:10])
            )
    raise FileNotFoundError(
        f"Cannot resolve proteome for {sample.sample}. Tried {rendered}"
        + (f" and pattern {pattern}" if pattern else "")
    )


def fasta_aliases(record: SeqRecord) -> set[str]:
    aliases = {record.id, record.name, record.description.split()[0]}
    for token in record.description.replace("|", " ").split():
        aliases.add(token)
    return {x for x in aliases if x}


def index_fasta(path: Path) -> dict[str, SeqRecord]:
    index: dict[str, SeqRecord] = {}
    for record in SeqIO.parse(str(path), "fasta"):
        for alias in fasta_aliases(record):
            index.setdefault(alias, record)
    if not index:
        raise ValueError(f"No FASTA records found in {path}")
    return index


def select_best_hits(
    grouped: dict[tuple[str, str], list[Hit]],
    ambiguity_ratio_cutoff: float,
    allow_ambiguous: bool,
) -> tuple[dict[tuple[str, str], Hit], list[dict[str, object]]]:
    selected: dict[tuple[str, str], Hit] = {}
    diagnostics: list[dict[str, object]] = []
    for key, candidates in sorted(grouped.items()):
        ranked = sorted(candidates, key=lambda h: (-h.normalized_score, -h.score, h.evalue, h.protein_id))
        best = ranked[0]
        second = ranked[1] if len(ranked) > 1 else None
        ratio = None
        ambiguous = False
        if second is not None and best.normalized_score > 0:
            ratio = second.normalized_score / best.normalized_score
            ambiguous = math.isfinite(ratio) and ratio >= ambiguity_ratio_cutoff
        if allow_ambiguous or not ambiguous:
            selected[key] = best
        diagnostics.append(
            {
                "sample": key[0],
                "ko": key[1],
                "n_significant_hits": len(ranked),
                "best_protein_id": best.protein_id,
                "best_score": best.score,
                "best_threshold": best.threshold,
                "best_normalized_score": best.normalized_score,
                "second_normalized_score": second.normalized_score if second else "",
                "second_to_best_ratio": ratio if ratio is not None else "",
                "ambiguous": ambiguous,
                "selected": allow_ambiguous or not ambiguous,
            }
        )
    return selected, diagnostics


def write_tsv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_command(command: Sequence[str], stdout: Path | None = None) -> None:
    print("+", " ".join(str(x) for x in command), file=sys.stderr)
    if stdout is None:
        subprocess.run(command, check=True)
    else:
        stdout.parent.mkdir(parents=True, exist_ok=True)
        with stdout.open("w", encoding="utf-8") as handle:
            subprocess.run(command, check=True, stdout=handle)


def ensure_tools(tools: Iterable[str]) -> None:
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        raise RuntimeError("Missing required executables: " + ", ".join(missing))


def build_concatenation(
    marker_alignments: list[tuple[str, Path]],
    samples: list[str],
    output_fasta: Path,
    partition_file: Path,
) -> list[dict[str, object]]:
    concatenated = {sample: [] for sample in samples}
    partitions: list[dict[str, object]] = []
    start = 1
    for ko, alignment_path in marker_alignments:
        records = {record.id: str(record.seq) for record in SeqIO.parse(str(alignment_path), "fasta")}
        lengths = {len(sequence) for sequence in records.values()}
        if len(lengths) != 1:
            raise ValueError(f"Alignment {alignment_path} has inconsistent sequence lengths")
        length = next(iter(lengths))
        end = start + length - 1
        for sample in samples:
            concatenated[sample].append(records.get(sample, "-" * length))
        partitions.append({"ko": ko, "start": start, "end": end, "length": length})
        start = end + 1
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(
        [SeqRecord(Seq("".join(concatenated[sample])), id=sample, description="") for sample in samples],
        str(output_fasta),
        "fasta",
    )
    with partition_file.open("w", encoding="utf-8") as handle:
        handle.write("#nexus\nbegin sets;\n")
        for row in partitions:
            handle.write(f"  charset {row['ko']} = {row['start']}-{row['end']};\n")
        handle.write("end;\n")
    return partitions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, default=Path("config/samples.tsv"))
    parser.add_argument("--results-root", type=Path, default=Path("results/per_genome"))
    parser.add_argument("--protein-template", default="results/per_genome/{sample}/proteins/proteins.faa")
    parser.add_argument("--proteome-root", type=Path, default=Path("../FungaltAI/data/proteome"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/phylogeny/ribosomal_markers"))
    parser.add_argument("--marker-prevalence", type=float, default=0.90)
    parser.add_argument("--final-marker-coverage", type=float, default=0.80)
    parser.add_argument("--min-genome-marker-coverage", type=float, default=0.70)
    parser.add_argument("--ambiguity-ratio", type=float, default=0.95)
    parser.add_argument("--allow-ambiguous", action="store_true")
    parser.add_argument("--min-length-ratio", type=float, default=0.50)
    parser.add_argument("--max-length-ratio", type=float, default=1.50)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--iqtree-model", default="LG+F+G4")
    parser.add_argument("--ufboot", type=int, default=1000)
    parser.add_argument("--alrt", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--skip-tree", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0 < args.marker_prevalence <= 1:
        raise ValueError("--marker-prevalence must be in (0, 1]")
    if not 0 < args.final_marker_coverage <= 1:
        raise ValueError("--final-marker-coverage must be in (0, 1]")
    if not 0 < args.min_genome_marker_coverage <= 1:
        raise ValueError("--min-genome-marker-coverage must be in (0, 1]")
    if not 0 < args.ambiguity_ratio <= 1:
        raise ValueError("--ambiguity-ratio must be in (0, 1]")

    samples = read_samples(args.samples)
    sample_names = [x.sample for x in samples]
    n_samples = len(samples)
    output_dir = args.output_dir
    raw_dir = output_dir / "markers_raw"
    aligned_dir = output_dir / "markers_aligned"
    trimmed_dir = output_dir / "markers_trimmed"
    qc_dir = output_dir / "qc"
    for directory in (raw_dir, aligned_dir, trimmed_dir, qc_dir):
        directory.mkdir(parents=True, exist_ok=True)

    hits_by_key: dict[tuple[str, str], list[Hit]] = defaultdict(list)
    marker_presence: dict[str, set[str]] = defaultdict(set)
    definitions: dict[str, Counter[str]] = defaultdict(Counter)
    for sample in samples:
        hit_path = args.results_root / sample.sample / "kofamscan" / "ribosome_significant_hits.tsv"
        for hit in read_hits(hit_path, sample.sample):
            hits_by_key[(sample.sample, hit.ko)].append(hit)
            marker_presence[hit.ko].add(sample.sample)
            if hit.definition:
                definitions[hit.ko][hit.definition] += 1

    prevalence_rows = []
    eligible_markers = set()
    for ko in sorted(marker_presence):
        prevalence = len(marker_presence[ko]) / n_samples
        if prevalence >= args.marker_prevalence:
            eligible_markers.add(ko)
        prevalence_rows.append(
            {
                "ko": ko,
                "definition": definitions[ko].most_common(1)[0][0] if definitions[ko] else "",
                "n_genomes_present": len(marker_presence[ko]),
                "n_expected_genomes": n_samples,
                "prevalence": f"{prevalence:.6f}",
                "passes_marker_prevalence": prevalence >= args.marker_prevalence,
            }
        )
    write_tsv(qc_dir / "marker_prevalence.tsv", prevalence_rows, ["ko", "definition", "n_genomes_present", "n_expected_genomes", "prevalence", "passes_marker_prevalence"])
    if not eligible_markers:
        raise RuntimeError("No KO marker passes the requested prevalence threshold")

    filtered_groups = {key: value for key, value in hits_by_key.items() if key[1] in eligible_markers}
    selected_hits, ambiguity_rows = select_best_hits(filtered_groups, args.ambiguity_ratio, args.allow_ambiguous)
    write_tsv(
        qc_dir / "best_hit_ambiguity.tsv",
        ambiguity_rows,
        ["sample", "ko", "n_significant_hits", "best_protein_id", "best_score", "best_threshold", "best_normalized_score", "second_normalized_score", "second_to_best_ratio", "ambiguous", "selected"],
    )

    selected_sequences: dict[str, list[SelectedSequence]] = defaultdict(list)
    extraction_rows = []
    for sample in samples:
        keys = [key for key in selected_hits if key[0] == sample.sample]
        if not keys:
            continue
        proteome_path = resolve_proteome(sample, args.protein_template, args.proteome_root)
        index = index_fasta(proteome_path)
        for key in sorted(keys):
            hit = selected_hits[key]
            record = index.get(hit.protein_id)
            extraction_rows.append({"sample": sample.sample, "ko": hit.ko, "protein_id": hit.protein_id, "proteome": str(proteome_path), "status": "found" if record is not None else "protein_id_not_found"})
            if record is None:
                continue
            ranked = sorted(filtered_groups[key], key=lambda h: (-h.normalized_score, -h.score, h.evalue, h.protein_id))
            second = ranked[1].normalized_score if len(ranked) > 1 else None
            ratio = second / hit.normalized_score if second is not None and hit.normalized_score > 0 else None
            sequence = str(record.seq).replace("*", "").upper()
            if sequence:
                selected_sequences[hit.ko].append(
                    SelectedSequence(
                        sample=sample.sample,
                        ko=hit.ko,
                        protein_id=hit.protein_id,
                        definition=hit.definition,
                        score=hit.score,
                        threshold=hit.threshold,
                        normalized_score=hit.normalized_score,
                        second_normalized_score=second,
                        ambiguity_ratio=ratio,
                        sequence=sequence,
                    )
                )
    write_tsv(qc_dir / "protein_extraction.tsv", extraction_rows, ["sample", "ko", "protein_id", "proteome", "status"])

    marker_qc_rows = []
    retained_markers: list[str] = []
    for ko in sorted(eligible_markers):
        records = selected_sequences.get(ko, [])
        lengths = [len(x.sequence) for x in records]
        median_length = statistics.median(lengths) if lengths else math.nan
        kept = [x for x in records if math.isfinite(median_length) and args.min_length_ratio * median_length <= len(x.sequence) <= args.max_length_ratio * median_length]
        coverage = len(kept) / n_samples
        retain = coverage >= args.final_marker_coverage and len(kept) >= 3
        if retain:
            retained_markers.append(ko)
            SeqIO.write(
                [SeqRecord(Seq(x.sequence), id=x.sample, description=f"ko={ko} protein_id={x.protein_id}") for x in kept],
                str(raw_dir / f"{ko}.faa"),
                "fasta",
            )
        marker_qc_rows.append(
            {
                "ko": ko,
                "definition": definitions[ko].most_common(1)[0][0] if definitions[ko] else "",
                "n_after_best_hit_filter": len(records),
                "median_length_aa": median_length,
                "n_after_length_filter": len(kept),
                "coverage_after_filters": f"{coverage:.6f}",
                "retained": retain,
            }
        )
    write_tsv(qc_dir / "marker_selection.tsv", marker_qc_rows, ["ko", "definition", "n_after_best_hit_filter", "median_length_aa", "n_after_length_filter", "coverage_after_filters", "retained"])
    (output_dir / "retained_markers.txt").write_text("\n".join(retained_markers) + "\n", encoding="utf-8")
    if len(retained_markers) < 5:
        raise RuntimeError(f"Only {len(retained_markers)} markers remain; inspect QC before tree inference")

    print(f"Samples={n_samples}; prevalent_KOs={len(eligible_markers)}; retained_markers={len(retained_markers)}", file=sys.stderr)
    if args.prepare_only:
        return

    ensure_tools(["mafft", "trimal"] + ([] if args.skip_tree else ["iqtree2"]))
    marker_alignments: list[tuple[str, Path]] = []
    alignment_rows = []
    for ko in retained_markers:
        raw = raw_dir / f"{ko}.faa"
        aligned = aligned_dir / f"{ko}.aln.faa"
        trimmed = trimmed_dir / f"{ko}.trimmed.faa"
        run_command(["mafft", "--auto", "--thread", str(args.threads), str(raw)], stdout=aligned)
        run_command(["trimal", "-in", str(aligned), "-out", str(trimmed), "-automated1"])
        records = list(SeqIO.parse(str(trimmed), "fasta"))
        lengths = {len(record.seq) for record in records}
        good = len(records) >= 3 and len(lengths) == 1 and next(iter(lengths), 0) > 0
        alignment_rows.append({"ko": ko, "retained_after_trimming": good, "n_sequences": len(records), "alignment_length": next(iter(lengths), 0) if len(lengths) == 1 else ""})
        if good:
            marker_alignments.append((ko, trimmed))
    write_tsv(qc_dir / "alignment_qc.tsv", alignment_rows, ["ko", "retained_after_trimming", "n_sequences", "alignment_length"])
    if len(marker_alignments) < 5:
        raise RuntimeError("Fewer than five valid trimmed marker alignments remain")

    alignment_presence = {
        ko: {record.id for record in SeqIO.parse(str(path), "fasta")}
        for ko, path in marker_alignments
    }
    genome_coverage_rows = []
    retained_sample_names = []
    for sample_name in sample_names:
        present = sum(sample_name in alignment_presence[ko] for ko, _ in marker_alignments)
        coverage = present / len(marker_alignments)
        retain = coverage >= args.min_genome_marker_coverage
        if retain:
            retained_sample_names.append(sample_name)
        genome_coverage_rows.append(
            {
                "sample": sample_name,
                "n_markers_present": present,
                "n_markers_total": len(marker_alignments),
                "marker_coverage": f"{coverage:.6f}",
                "retained": retain,
            }
        )
    write_tsv(
        qc_dir / "genome_marker_coverage.tsv",
        genome_coverage_rows,
        ["sample", "n_markers_present", "n_markers_total", "marker_coverage", "retained"],
    )
    if len(retained_sample_names) < 3:
        raise RuntimeError("Fewer than three genomes pass the marker-coverage threshold")

    concatenated = output_dir / "ribosomal_markers_concatenated.faa"
    partitions_path = output_dir / "ribosomal_markers_partitions.nex"
    partition_rows = build_concatenation(
        marker_alignments,
        retained_sample_names,
        concatenated,
        partitions_path,
    )
    write_tsv(qc_dir / "partition_coordinates.tsv", partition_rows, ["ko", "start", "end", "length"])
    if args.skip_tree:
        return

    prefix = output_dir / "fungal_ribosomal_species_tree"
    command = [
        "iqtree2", "-s", str(concatenated), "-p", str(partitions_path),
        "-m", args.iqtree_model, "-T", str(args.threads), "--seed", str(args.seed),
        "--prefix", str(prefix),
    ]
    if args.ufboot >= 1000:
        command.extend(["-B", str(args.ufboot)])
    if args.alrt >= 1000:
        command.extend(["--alrt", str(args.alrt)])
    run_command(command)
    print(f"Final tree: {prefix}.treefile", file=sys.stderr)


if __name__ == "__main__":
    main()
