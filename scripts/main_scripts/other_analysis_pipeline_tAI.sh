#!/bin/bash

#This script is the second and final part of the process for calculating the tAI
# (transfer RNA Adaptation Index) for each sequence in the CDS (coding sequence) folder.
#In the previous script, we performed the tRNAscan-SE computation.
# The overall process has been divided into two smaller scripts to optimize performance,
# as tRNAscan-SE can take a significant amount of time to run.

#First part of script is: /home/norbert_s/tAI/all_chunks/tRNAscanse_pipeline_tAI.sh

# ================================================
# Script to Calculate tRNA Adaptation Index (tAI)
#
# This bash script automates the process of calculating the
# tRNA Adaptation Index (tAI) for coding sequences (CDS)
# in a genomic dataset. It integrates tools, like:
# tRNAscan-SE, CodonM, and CodonW.
#
# Key Steps:
# 1. **File Setup**:
#    - Creates temporary folders and checks if the necessary
#      input directories (genomic_fna, cds) exist.
#    - Builds a mapping of genomic sequence files using common
#      prefixes (assembly IDs).
#
# 2. **Data Pairing**:
#    - Pairs genomic sequences with their corresponding CDS
#      sequences based on shared prefixes.
#
# 3. **Processing Loop**:
#    - For each directory (chunk), the script:
#      - Runs tRNAscan-SE to identify tRNA genes.
#      - Executes CodonM for codon usage analysis.
#      - Runs CodonW for protein sequence statistics.
#      - Merges first part of fasta header, because CodonW process
#       only first 20 signs from fasta header, what can make output non specific
#      - Merges protein IDs with CodonW output.
#      - Computes tAI using an R script based on the gathered data.
#
# 4. **Output**:
#    - Generates per-protein tAI values and a summary file for
#      each dataset chunk.
#
# Dependencies:
# - Python (`count_genes_coding_tRNA.py`)
# - Perl (`CodonM`)
# - R (`tAI.R`)
# - tRNAscan-SE (for tRNA prediction)
#
# Usage:
# To run the script, execute the following command:
# bash tAI_calculation_script.sh
# ================================================



# function to check for errors after each command
check_error() {
    if [ $? -ne 0 ]; then
        echo "Error in $1"
        exit 1
    fi
}


#setting temporary folder
mkdir -p /home/norbert_s/tmp
export TMPDIR=/home/norbert_s/tmp

#setting pathways to genomic and cds seqs
genomic_folder_path="/db/183_fungi_mai_2022/genomic_fna"
cds_folder_path="/home/norbert_s/tAI/all_chunks/cds"

# genomic_folder_path="/home/norbert_s/genome"
# cds_folder_path="/home/norbert_s/183_CDS"

#checking if input folders exist
if [ ! -d "$genomic_folder_path"  ]; then
    echo "$genomic_folder_path doesn't exist!"
    exit 1
fi

if [ ! -d "$cds_folder_path"  ]; then
    echo "$cds_folder_path doesn't exist!"
    exit 1
fi



# Creating arrays with genomic and CDS sequence files
genomic_seq_files=($(find "$genomic_folder_path" -type f -name "*.fna"))
cds_seq_files=($(find "$cds_folder_path" -type f -name "*.fna"))

tRNAscan_SE_files=($(find "$tRNAscan_SE_output" -type f -name "*.fna"))

#creating table to store files with prefix
declare -A genomic_files_map

#initialising map 
echo "Creating map started"
for genomic_seq in "${genomic_seq_files[@]}"; do
    genomic_prefix=$(basename "$genomic_seq" | cut -c1-15) #creating prefix whis is common with assemblyID
    genomic_files_map["$genomic_prefix"]+="$genomic_seq"
done
echo "Creating map done"    


#creating pairs with common prefix between cds folder and
file_pairs=()
echo "Creating pairs started"
for cds_seq in "${cds_seq_files[@]}"; do
    cds_prefix=$(basename "$cds_seq" | cut -c1-15) #creating prefix whis is common with assemblyID
        if [[ -n "${genomic_files_map[$cds_prefix]}" ]]; then
        for genomic_seq in ${genomic_files_map[$cds_prefix]}; do #combining genomic seq and cds seq to pairs according to common prefix (assemblyID)
            file_pairs+=("$genomic_seq $cds_seq")
        done
    fi
done
echo "Creating pairs done"


base_folder="/home/norbert_s/tAI/all_chunks/chunk"
script_path="/home/norbert_s/tAI"

#iterating through ten folders where are written datas from previous step
for folder_index in {1..10}; do
    echo "$folder_index"
    tRNAscan_SE_folder="${base_folder}${folder_index}/tRNAscan_SE"
    counted_genes_coding_tRNA_folder="${base_folder}${folder_index}/counted_genes_coding_tRNA_folder"   
    codonm_folder="${base_folder}${folder_index}/codonm_folder"
    codonw_folder="${base_folder}${folder_index}/codonw_folder"
    tAI_folder="${base_folder}${folder_index}/tAI_folder"

    mkdir -p "$counted_genes_coding_tRNA_folder"
    mkdir -p "$codonm_folder"
    mkdir -p "$codonw_folder"
    mkdir -p "$tAI_folder"

    if [ -d "$tRNAscan_SE_folder" ]; then
        for file_path in "$tRNAscan_SE_folder"/*; do
            if [ -f "$file_path" ]; then
                prefix=$(basename "$file_path" | tail -c 14)

                echo "Prefix $prefix"

                #Python script counting genes that code tRNA (my own script)
                echo ""
                echo "Python script count_genes_coding_tRNA.py started for $file_path"
                python3 "$script_path/count_genes_coding_tRNA.py" "$file_path" "$counted_genes_coding_tRNA_folder/counted_genes_coding_tRNA_$prefix.txt"
                check_error "count_genes_coding_tRNA.py for $prefix"

                #CodonM is Mardo Dos Reis's script (Perl) that generates a matrix of codon frequencies per ORF
                echo " "
                echo "CodonM started for $prefix"
                perl "$script_path/codonM" "$(ls $cds_folder_path/$prefix* | head -n 1)" "$codonm_folder/codonM_output_$prefix"
                check_error "codonM for $prefix"

                #CodonW counts some usefull statistics in all fasta (protein) sequnces.
                echo " "
                echo "CodonW started for $prefix"
                codonw "$(find "$cds_folder_path" -name "$prefix*" -print -quit)" "$codonw_folder/codonw_$prefix" "$codonw_folder/bulk_$prefix" -nomenu -nowarn -silent -enc -gc -gc3s -L_aa -noblk
                check_error "Codonw for $prefix"
#				sed -i -E 's/^[[:space:]]+/ /; s/[[:space:]]+/        /g' "$codonw_folder/codonw_$prefix" #replace multi-spaces signs by one tab (not user-friendly codonw)
				head "$codonw_folder/codonw_$prefix"

#                protein_id_column=$(find "$cds_folder_path" -name "$prefix*" -print | head -n 1)  # Get one file matching the prefix

#                if [ -n "$protein_id_column" ]; then
#                    # Assuming protein_id is in the form of [protein_id=...], extract them
#                    protein_ids=$(grep -o '\[protein_id=[^]]*\]' "$protein_id_column" | sed 's/\[protein_id=//; s/\]//')
#                    echo "${protein_ids[@]}"
#                    # Add header for Protein_ID
#                    protein_ids="Protein_id"$'\n'"$protein_ids"
#
#                    # Get the number of lines in both files to compare
#                    codonw_line_count=$(wc -l < "$codonw_folder/codonw_$prefix")
#                    protein_id_line_count=$(echo "$protein_ids" | wc -l)
#
#                    # Compare the number of lines to ensure they match
#                    if [ "$codonw_line_count" -eq "$protein_id_line_count" ]; then
#                        echo "Line counts match. Proceeding with merging."
#
#                        # Create a temporary file to store the merged content
#                        tmp_file=$(mktemp)
#
#                        paste <(echo "${protein_ids[@]}") "$codonw_folder/codonw_$prefix" > "$tmp_file"
#
#
##                        awk '{print $1 "\t" $2 "\t" $7 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' "$tmp_file" > file1.txt && mv file1.txt "$tmp_file"
#
#
#                        # Move the temporary file to the final destination
#                        mv "$tmp_file" "$codonw_folder/codonw_$prefix"
#                    else
#                        echo "Warning: Line counts do not match. Cannot merge files."
#                    fi
#                else
#                    echo "Error: No files found matching prefix $prefix."
#                fi

                #Script basis on the Mario Dos Reis https://github.com/mariodosreis/tai counting tAI
                echo "Rscript started for $prefix"
                echo "output folder: $tAI_folder"
                Rscript "$script_path/tAI.R" "$counted_genes_coding_tRNA_folder/counted_genes_coding_tRNA_$prefix.txt" "$codonm_folder/codonM_output_$prefix" "$codonw_folder/codonw_$prefix" $tAI_folder/each_protein_tAI_$prefix.txt $tAI_folder/summary_tAI_$prefix.txt
                check_error "Rscript tAI.R for $prefix in folder number $folder_index"



                echo " "
                echo " "

			fi
        done

    fi
done