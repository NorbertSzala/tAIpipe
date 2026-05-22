#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage:
  ./create_test_dataset.sh -C path/to/all/CDS -G path/to/all/genomes -A "GCA_000219625.1,GCA_000006255.1"

Required:
  -C    Path to directory with all CDS FASTA files
  -G    Path to directory with all genome FASTA files
  -A    Accessions to downsample, comma-separated or space-separated

Optional:
  -O    Output test data directory [default: resources/test_data]
  -L    Max bp from first genome contig [default: 1000000]
  -N    Fallback max CDS records if contig matching fails [default: 1000]

Output:
  resources/test_data/
  ├── CDS/
  │   └── ACCESSION_first_contig_CDS.fna
  └── genome/
      └── ACCESSION_first_contig_1Mbp.fna

Example:
  ./create_test_dataset.sh \\
    -C data/all_CDS \\
    -G data/all_genomes \\
    -A "GCA_000219625.1,GCA_000006255.1"
EOF
}

# ------------------------
# --- Argument Parsing ---
# ------------------------

CDS_DIR=""
GENOME_DIR=""
ACCESSIONS_RAW=""
OUTDIR="resources/test_data"
MAX_BP=1000000
MAX_CDS=1000

while getopts ":C:G:A:O:L:N:h" opt; do
    case "${opt}" in
        C) CDS_DIR="${OPTARG}" ;;
        G) GENOME_DIR="${OPTARG}" ;;
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

if [[ -z "${CDS_DIR}" || -z "${GENOME_DIR}" || -z "${ACCESSIONS_RAW}" ]]; then
    echo "ERROR: Missing required arguments." >&2
    usage
    exit 1
fi

if [[ ! -d "${CDS_DIR}" ]]; then
    echo "ERROR: CDS directory does not exist: ${CDS_DIR}" >&2
    exit 1
fi

if [[ ! -d "${GENOME_DIR}" ]]; then
    echo "ERROR: Genome directory does not exist: ${GENOME_DIR}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}/CDS" "${OUTDIR}/genome"

ACCESSIONS_RAW="${ACCESSIONS_RAW//,/ }"
read -r -a ACCESSIONS <<< "${ACCESSIONS_RAW}"


# ------------------
# --- Find files ---
# ------------------

find_single_file() {
    local dir="$1"
    local accession="$2"

    mapfile -t matches < <(
        find "${dir}" -maxdepth 1 -type f \
            \( \
                -name "${accession}*.fna" -o \
                -name "${accession}*.fa" -o \
                -name "${accession}*.fasta" -o \
                -name "${accession}*.fna.gz" -o \
                -name "${accession}*.fa.gz" -o \
                -name "${accession}*.fasta.gz" \
            \) \
            | sort
    )

    if [[ "${#matches[@]}" -eq 0 ]]; then
        echo "ERROR: No FASTA file found for accession ${accession} in ${dir}" >&2
        return 1
    fi

    if [[ "${#matches[@]}" -gt 1 ]]; then
        echo "ERROR: Multiple FASTA files found for accession ${accession} in ${dir}:" >&2
        printf '%s\n' "${matches[@]}" >&2
        return 1
    fi

    echo "${matches[0]}"
}


downsample_one() {
    local accession="$1"
    local genome_in="$2"
    local cds_in="$3"
    local genome_out="$4"
    local cds_out="$5"

    python3 - "$accession" "$genome_in" "$cds_in" "$genome_out" "$cds_out" "$MAX_BP" "$MAX_CDS" <<'PY'
import gzip
import sys
from pathlib import Path

accession = sys.argv[1]
genome_in = Path(sys.argv[2])
cds_in = Path(sys.argv[3])
genome_out = Path(sys.argv[4])
cds_out = Path(sys.argv[5])
max_bp = int(sys.argv[6])
max_cds = int(sys.argv[7])


def open_maybe_gzip(path, mode="rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def fasta_iter(path):
    header = None
    seq_chunks = []

    with open_maybe_gzip(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)

                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

    if header is not None:
        yield header, "".join(seq_chunks)


def wrap(seq, width=80):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def write_fasta_record(handle, header, seq):
    handle.write(f">{header}\n")
    for chunk in wrap(seq):
        handle.write(chunk + "\n")


# Extract first genome contig and truncate to max_bp
try:
    first_header, first_seq = next(fasta_iter(genome_in))
except StopIteration:
    raise SystemExit(f"ERROR: empty genome FASTA: {genome_in}")

first_contig_id = first_header.split()[0]
first_seq_short = first_seq[:max_bp]

genome_out.parent.mkdir(parents=True, exist_ok=True)
cds_out.parent.mkdir(parents=True, exist_ok=True)

with open(genome_out, "w") as out:
    write_fasta_record(
        out,
        f"{first_contig_id}|{accession}|first_contig_first_{len(first_seq_short)}bp",
        first_seq_short,
    )


# Extract CDS corresponding to first contig if CDS headers contain contig ID.
# If no such records are found, fallback to first max_cds CDS records.
matched = []
fallback = []

for h, s in fasta_iter(cds_in):
    if len(fallback) < max_cds:
        fallback.append((h, s))

    if first_contig_id in h:
        matched.append((h, s))

if matched:
    selected = matched
    mode = "contig_header_match"
else:
    selected = fallback
    mode = f"fallback_first_{max_cds}_cds"

with open(cds_out, "w") as out:
    for h, s in selected:
        write_fasta_record(out, h, s)


print(f"accession: {accession}")
print(f"first_contig_id: {first_contig_id}")
print(f"genome_input: {genome_in}")
print(f"cds_input: {cds_in}")
print(f"genome_output: {genome_out}")
print(f"cds_output: {cds_out}")
print(f"genome_bp_written: {len(first_seq_short)}")
print(f"cds_records_written: {len(selected)}")
print(f"cds_selection_mode: {mode}")
print()
print("-----------------------------------------------")
print()
PY
}


for accession in "${ACCESSIONS[@]}"; do
    echo "Processing accession: ${accession}" >&2

    genome_in="$(find_single_file "${GENOME_DIR}" "${accession}")"
    cds_in="$(find_single_file "${CDS_DIR}" "${accession}")"

    genome_out="${OUTDIR}/genome/${accession}_first_contig_1Mbp.fna"
    cds_out="${OUTDIR}/CDS/${accession}_first_contig_CDS.fna"

    downsample_one \
        "${accession}" \
        "${genome_in}" \
        "${cds_in}" \
        "${genome_out}" \
        "${cds_out}"
done

echo ""
echo "Done."
echo "Test genomes: ${OUTDIR}/genome"
echo "Test CDS:     ${OUTDIR}/CDS"
echo "--------------------------------"