#!/bin/bash
# =====================================================================
#  Enhanced tRNAscan-SE pipeline — version 2.0
# =====================================================================
#  Author: Norbert S.
#  Purpose:
#    Automates tRNAscan-SE analysis for genomic sequences grouped
#    into chunks, matching CDS files by assembly ID prefix.
#
#  Features:
#    ✓ Automatic pairing of genome and CDS files by prefix
#    ✓ Parallel execution with CPU usage control
#    ✓ Organized output folders (10 chunks)
#    ✓ Error and progress logging
#    ✓ Resilient to partial failures or duplicated IDs
#
#  Second stage: /home/norbert_s/tAI/all_chunks/other_analysis_pipeline_tAI.sh
# =====================================================================

set -euo pipefail
IFS=$'\n\t'

# =========================== CONFIGURATION ===========================
genomic_folder_path="/db/183_fungi_mai_2022/genomic_fna"
cds_folder_path="/home/norbert_s/tAI/all_chunks/cds"
output_base_folder="output_shorted_pipeline_tAI"
files_per_folder=19
max_jobs=6  # Limit of parallel jobs
export TMPDIR="/home/norbert_s/tmp"
log_file="pipeline_log.txt"
# =====================================================================

# --------------------- Helper Functions ---------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$log_file"
}

check_error() {
    if [ $? -ne 0 ]; then
        log "Error in $1"
        exit 1
    fi
}

# Ensure required folders exist
if [ ! -d "$genomic_folder_path" ]; then
    log "Folder not found: $genomic_folder_path"
    exit 1
fi

if [ ! -d "$cds_folder_path" ]; then
    log "Folder not found: $cds_folder_path"
    exit 1
fi

# Create output folders
for i in {1..10}; do
    mkdir -p "${output_base_folder}${i}"
done
log "Output folders created."

# --------------------- Pair Genomic & CDS ---------------------
log "Scanning genomic and CDS files..."

mapfile -t genomic_seq_files < <(find "$genomic_folder_path" -type f -name "*.fna")
mapfile -t cds_seq_files < <(find "$cds_folder_path" -type f -name "*.fna")

declare -A genomic_files_map
declare -A seen_prefixes

# Map genomic files by assembly prefix
for genomic_seq in "${genomic_seq_files[@]}"; do
    prefix=$(basename "$genomic_seq" | cut -c1-13)
    genomic_files_map["$prefix"]="$genomic_seq"
done

file_pairs=()
for cds_seq in "${cds_seq_files[@]}"; do
    prefix=$(basename "$cds_seq" | cut -c1-13)
    if [[ -n "${genomic_files_map[$prefix]:-}" ]]; then
        # Avoid duplicates
        if [[ -z "${seen_prefixes[$prefix]:-}" ]]; then
            file_pairs+=("${genomic_files_map[$prefix]}::${cds_seq}")
            seen_prefixes["$prefix"]=1
        fi
    fi
done

num_pairs=${#file_pairs[@]}
log "Created $num_pairs valid genome/CDS pairs."

# --------------------- Process Each Pair ---------------------
process_pair() {
    pair=$1
    counter=$2
    output_folder=$3

    genomic_seq=${pair%%::*}
    cds_seq=${pair##*::}
    prefix=$(basename "$genomic_seq" | cut -c1-13)

    log "[${counter}] Running tRNAscan-SE for prefix: $prefix"

    output_file="$output_folder/tRNAscan-SE_output_$prefix"
    tRNAscan-SE -D -o "$output_file" "$genomic_seq" > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        log "tRNAscan-SE failed for $prefix"
    else
        log "Completed $prefix → saved to $output_file"
    fi
}
export -f process_pair
export -f log
export output_base_folder

# --------------------- Parallel Execution ---------------------
log "Starting processing..."
start_time=$(date +%s)

for ((i=0; i<num_pairs; i++)); do
    folder_index=$(( (i / files_per_folder) + 1 ))
    output_folder="${output_base_folder}${folder_index}"

    process_pair "${file_pairs[$i]}" $((i+1)) "$output_folder" &

    # Limit number of concurrent jobs
    while (( $(jobs -r | wc -l) >= max_jobs )); do
        sleep 2
    done
done

wait
end_time=$(date +%s)
runtime=$((end_time - start_time))
log "All $num_pairs pairs processed in ${runtime}s."
