#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Description:
    Create a consistent test subset containing genome, CDS, proteome and GFF
    files for selected accessions.

    The first genome contig is truncated to a requested length. CDS and GFF
    records are retained only when they belong to this contig and fit within
    the retained interval. The proteome is filtered using protein_id values
    extracted from the selected CDS records.

Usage:
    create_test_dataset.sh \
        -C CDS_DIR \
        -G GENOME_DIR \
        -P PROTEOME_DIR \
        -F GFF_DIR \
        -A ACCESSIONS \
        [-O OUTPUT_DIR] \
        [-L MAX_BP] \
        [-N MAX_CDS]

Required:
    -C    Directory containing CDS FASTA files
    -G    Directory containing genome FASTA files
    -P    Directory containing protein FASTA files
    -F    Directory containing GFF/GFF3 files
    -A    Comma- or space-separated accessions

Optional:
    -O    Output directory [default: tests/data]
    -L    Maximum retained length of first contig [default: 1000000]
    -N    Maximum fallback number of CDS records [default: 1000]

Output:
    tests/data/
    ├── CDS/
    ├── genome/
    ├── proteome/
    └── gff/
EOF
}

CDS_DIR=""
GENOME_DIR=""
PROTEOME_DIR=""
GFF_DIR=""
ACCESSIONS_RAW=""
OUTDIR="tests/data"
MAX_BP=1000000
MAX_CDS=1000

while getopts ":C:G:P:F:A:O:L:N:h" opt; do
    case "${opt}" in
        C) CDS_DIR="${OPTARG}" ;;
        G) GENOME_DIR="${OPTARG}" ;;
        P) PROTEOME_DIR="${OPTARG}" ;;
        F) GFF_DIR="${OPTARG}" ;;
        A) ACCESSIONS_RAW="${OPTARG}" ;;
        O) OUTDIR="${OPTARG}" ;;
        L) MAX_BP="${OPTARG}" ;;
        N) MAX_CDS="${OPTARG}" ;;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "ERROR: Invalid option: -${OPTARG}" >&2
            usage
            exit 1
            ;;
        :)
            echo "ERROR: Option -${OPTARG} requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

for value_name in \
    CDS_DIR GENOME_DIR PROTEOME_DIR GFF_DIR ACCESSIONS_RAW
do
    if [[ -z "${!value_name}" ]]; then
        echo "ERROR: Missing required argument: ${value_name}" >&2
        usage
        exit 1
    fi
done

for directory in \
    "${CDS_DIR}" \
    "${GENOME_DIR}" \
    "${PROTEOME_DIR}" \
    "${GFF_DIR}"
do
    if [[ ! -d "${directory}" ]]; then
        echo "ERROR: Directory does not exist: ${directory}" >&2
        exit 1
    fi
done

mkdir -p \
    "${OUTDIR}/CDS" \
    "${OUTDIR}/genome" \
    "${OUTDIR}/proteome" \
    "${OUTDIR}/gff"

ACCESSIONS_RAW="${ACCESSIONS_RAW//,/ }"
read -r -a ACCESSIONS <<< "${ACCESSIONS_RAW}"


find_single_fasta() {
    local directory="$1"
    local accession="$2"
    local data_type="$3"

    local patterns=()

    case "${data_type}" in
        genome|cds)
            patterns=(
                "${accession}*.fna"
                "${accession}*.fa"
                "${accession}*.fasta"
                "${accession}*.fna.gz"
                "${accession}*.fa.gz"
                "${accession}*.fasta.gz"
            )
            ;;
        protein)
            patterns=(
                "${accession}*.faa"
                "${accession}*.pep"
                "${accession}*.fa"
                "${accession}*.fasta"
                "${accession}*.faa.gz"
                "${accession}*.pep.gz"
                "${accession}*.fa.gz"
                "${accession}*.fasta.gz"
            )
            ;;
        *)
            echo "ERROR: Unknown FASTA type: ${data_type}" >&2
            return 1
            ;;
    esac

    local find_args=()

    for pattern in "${patterns[@]}"; do
        if [[ "${#find_args[@]}" -gt 0 ]]; then
            find_args+=("-o")
        fi
        find_args+=("-name" "${pattern}")
    done

    mapfile -t matches < <(
        find "${directory}" -maxdepth 1 -type f \
            \( "${find_args[@]}" \) |
            sort
    )

    if [[ "${#matches[@]}" -eq 0 ]]; then
        echo "ERROR: No ${data_type} FASTA for ${accession}" >&2
        return 1
    fi

    if [[ "${#matches[@]}" -gt 1 ]]; then
        echo "ERROR: Multiple ${data_type} FASTA files for ${accession}:" >&2
        printf '%s\n' "${matches[@]}" >&2
        return 1
    fi

    printf '%s\n' "${matches[0]}"
}


find_single_gff() {
    local directory="$1"
    local accession="$2"

    mapfile -t matches < <(
        find "${directory}" -maxdepth 1 -type f \
            \( \
                -name "${accession}*.gff" -o \
                -name "${accession}*.gff3" -o \
                -name "${accession}*.gff.gz" -o \
                -name "${accession}*.gff3.gz" \
            \) |
            sort
    )

    if [[ "${#matches[@]}" -eq 0 ]]; then
        echo "ERROR: No GFF file for ${accession}" >&2
        return 1
    fi

    if [[ "${#matches[@]}" -gt 1 ]]; then
        echo "ERROR: Multiple GFF files for ${accession}:" >&2
        printf '%s\n' "${matches[@]}" >&2
        return 1
    fi

    printf '%s\n' "${matches[0]}"
}


downsample_one() {
    local accession="$1"
    local genome_in="$2"
    local cds_in="$3"
    local proteome_in="$4"
    local gff_in="$5"
    local genome_out="$6"
    local cds_out="$7"
    local proteome_out="$8"
    local gff_out="$9"

    python3 - \
        "${accession}" \
        "${genome_in}" \
        "${cds_in}" \
        "${proteome_in}" \
        "${gff_in}" \
        "${genome_out}" \
        "${cds_out}" \
        "${proteome_out}" \
        "${gff_out}" \
        "${MAX_BP}" \
        "${MAX_CDS}" <<'PY'
import gzip
import re
import sys
from pathlib import Path


(
    accession,
    genome_in,
    cds_in,
    proteome_in,
    gff_in,
    genome_out,
    cds_out,
    proteome_out,
    gff_out,
    max_bp,
    max_cds,
) = sys.argv[1:]

genome_in = Path(genome_in)
cds_in = Path(cds_in)
proteome_in = Path(proteome_in)
gff_in = Path(gff_in)

genome_out = Path(genome_out)
cds_out = Path(cds_out)
proteome_out = Path(proteome_out)
gff_out = Path(gff_out)

max_bp = int(max_bp)
max_cds = int(max_cds)

PROTEIN_ID_PATTERN = re.compile(r"\[protein_id=([^\]]+)\]")
LOCATION_PATTERN = re.compile(r"\[location=([^\]]+)\]")


def open_maybe_gzip(path, mode="rt"):
    path = Path(path)

    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")

    return path.open(mode, encoding="utf-8")


def fasta_iter(path):
    header = None
    sequence_parts = []

    with open_maybe_gzip(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)

                header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())

    if header is not None:
        yield header, "".join(sequence_parts)


def write_fasta_record(handle, header, sequence, width=80):
    handle.write(f">{header}\n")

    for start in range(0, len(sequence), width):
        handle.write(sequence[start:start + width] + "\n")


def normalize_contig_id(value):
    value = value.strip()

    if value.startswith("lcl|"):
        value = value[4:]

    return value


def cds_contig_id(header):
    token = header.split()[0]
    token = normalize_contig_id(token)

    if "_cds_" in token:
        token = token.split("_cds_", maxsplit=1)[0]

    return token


def cds_fits_interval(header, retained_length):
    location_match = LOCATION_PATTERN.search(header)

    if location_match is None:
        return True

    coordinates = [
        int(value)
        for value in re.findall(r"\d+", location_match.group(1))
    ]

    if not coordinates:
        return True

    return min(coordinates) >= 1 and max(coordinates) <= retained_length


def extract_protein_id(header):
    match = PROTEIN_ID_PATTERN.search(header)
    return match.group(1).strip() if match else None


try:
    first_genome_header, first_genome_sequence = next(
        fasta_iter(genome_in)
    )
except StopIteration:
    raise SystemExit(f"ERROR: Empty genome FASTA: {genome_in}")

first_contig_id = first_genome_header.split()[0]
normalized_first_contig = normalize_contig_id(first_contig_id)

retained_sequence = first_genome_sequence[:max_bp]
retained_length = len(retained_sequence)

for output_path in (
    genome_out,
    cds_out,
    proteome_out,
    gff_out,
):
    output_path.parent.mkdir(parents=True, exist_ok=True)

with genome_out.open("w", encoding="utf-8") as handle:
    write_fasta_record(
        handle,
        (
            f"{first_contig_id}|{accession}|"
            f"first_{retained_length}bp"
        ),
        retained_sequence,
    )


matched_cds = []
fallback_cds = []

for header, sequence in fasta_iter(cds_in):
    if len(fallback_cds) < max_cds:
        fallback_cds.append((header, sequence))

    same_contig = (
        normalize_contig_id(cds_contig_id(header))
        == normalized_first_contig
    )

    if same_contig and cds_fits_interval(
        header,
        retained_length,
    ):
        matched_cds.append((header, sequence))

if matched_cds:
    selected_cds = matched_cds
    selection_mode = "first_contig_within_retained_interval"
else:
    selected_cds = fallback_cds
    selection_mode = f"fallback_first_{len(selected_cds)}_cds"
    print(
        "WARNING: No CDS records matched the retained genome interval; "
        "using fallback CDS records.",
        file=sys.stderr,
    )

selected_protein_ids = set()

with cds_out.open("w", encoding="utf-8") as handle:
    for header, sequence in selected_cds:
        write_fasta_record(handle, header, sequence)

        protein_id = extract_protein_id(header)

        if protein_id:
            selected_protein_ids.add(protein_id)

if not selected_protein_ids:
    raise SystemExit(
        "ERROR: No protein_id attributes found in selected CDS records"
    )


selected_proteins = []
observed_protein_ids = set()

for header, sequence in fasta_iter(proteome_in):
    protein_id = header.split()[0]

    if protein_id in selected_protein_ids:
        selected_proteins.append((header, sequence))
        observed_protein_ids.add(protein_id)

missing_proteins = selected_protein_ids - observed_protein_ids

if missing_proteins:
    preview = ", ".join(sorted(missing_proteins)[:20])

    raise SystemExit(
        f"ERROR: {len(missing_proteins)} selected CDS proteins are "
        f"absent from the proteome. Examples: {preview}"
    )

with proteome_out.open("w", encoding="utf-8") as handle:
    for header, sequence in selected_proteins:
        write_fasta_record(handle, header, sequence)


gff_records_written = 0

with (
    open_maybe_gzip(gff_in) as source,
    gff_out.open("w", encoding="utf-8") as target,
):
    target.write("##gff-version 3\n")
    target.write(
        f"##sequence-region {first_contig_id} "
        f"1 {retained_length}\n"
    )

    for raw_line in source:
        if raw_line.startswith("##FASTA"):
            break

        if raw_line.startswith("#"):
            continue

        fields = raw_line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        seqid = normalize_contig_id(fields[0])

        try:
            start = int(fields[3])
            end = int(fields[4])
        except ValueError:
            continue

        if (
            seqid == normalized_first_contig
            and start >= 1
            and end <= retained_length
        ):
            target.write(raw_line)
            gff_records_written += 1

if gff_records_written == 0:
    print(
        "WARNING: No GFF records matched the retained interval.",
        file=sys.stderr,
    )


print(f"accession: {accession}")
print(f"first_contig_id: {first_contig_id}")
print(f"retained_genome_bp: {retained_length}")
print(f"cds_selection_mode: {selection_mode}")
print(f"selected_cds_records: {len(selected_cds)}")
print(f"selected_protein_records: {len(selected_proteins)}")
print(f"selected_gff_records: {gff_records_written}")
print(f"genome_output: {genome_out}")
print(f"cds_output: {cds_out}")
print(f"proteome_output: {proteome_out}")
print(f"gff_output: {gff_out}")
print()
PY
}


for accession in "${ACCESSIONS[@]}"; do
    echo "Processing accession: ${accession}" >&2

    genome_in="$(
        find_single_fasta "${GENOME_DIR}" "${accession}" genome
    )"

    cds_in="$(
        find_single_fasta "${CDS_DIR}" "${accession}" cds
    )"

    proteome_in="$(
        find_single_fasta "${PROTEOME_DIR}" "${accession}" protein
    )"

    gff_in="$(
        find_single_gff "${GFF_DIR}" "${accession}"
    )"

    genome_out=(
        "${OUTDIR}/genome/"
        "${accession}_first_contig_1Mbp.fna"
    )
    genome_out="${genome_out[*]}"

    cds_out=(
        "${OUTDIR}/CDS/"
        "${accession}_first_contig_CDS.fna"
    )
    cds_out="${cds_out[*]}"

    proteome_out=(
        "${OUTDIR}/proteome/"
        "${accession}_first_contig_protein.faa"
    )
    proteome_out="${proteome_out[*]}"

    gff_out=(
        "${OUTDIR}/gff/"
        "${accession}_first_contig_1Mbp.gff3"
    )
    gff_out="${gff_out[*]}"

    downsample_one \
        "${accession}" \
        "${genome_in}" \
        "${cds_in}" \
        "${proteome_in}" \
        "${gff_in}" \
        "${genome_out}" \
        "${cds_out}" \
        "${proteome_out}" \
        "${gff_out}"
done

echo
echo "Done."
echo "Genome:  ${OUTDIR}/genome"
echo "CDS:     ${OUTDIR}/CDS"
echo "Proteome:${OUTDIR}/proteome"
echo "GFF:     ${OUTDIR}/gff"