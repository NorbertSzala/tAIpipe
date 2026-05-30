#!/usr/bin/env python

"""
Description:
    Aggregates per-gene metrics from cubar summary files with protein 
    structural features (length, calculated LCRs with 3 and 10 bins profiling)
    and phenotypic metadata to create the master dataset table.
    Dynamically loads genetic code overrides from tRNAscan-SE .gcode files.

Usage:
    python compile_metadata_samples.py -D data/tutorial_data/input/metadata/dataset_test.tsv \
        -P data/tutorial_data/output/per_genome -C data/tutorial_data/input/CDS \
        -G data/tutorial_data/input/trnascanse -O data/tutorial_data/input/metadata/samples_test.tsv
"""

from pathlib import Path
import argparse
import pandas as pd
import math
from collections import Counter

# -----------------
# --- Arguments ---
def arguments():
    p = argparse.ArgumentParser(
        prog="compile_metadata_samples",
        description="Compile master metadata table linking translation metrics with structural protein profiles",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("-D", "--dataset-sheet", help="Path to input dataset_test.tsv steering sheet", required=True, type=Path)
    p.add_argument("-P", "--per-genome-dir", help="Path to per_genome output directory containing cubar summaries", required=True, type=Path)
    p.add_argument("-C", "--cds-dir", help="Path to directory containing input CDS FASTA files", required=True, type=Path)
    p.add_argument("-G", "--gcode-dir", help="Path to directory containing tRNAscan-SE .gcode files", required=True, type=Path)
    p.add_argument("-O", "--output", help="Path to output master samples_test.tsv file", required=True, type=Path)
    return p.parse_args()

# ------------------------
# --- Helper functions ---
# ------------------------

STANDARD_GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def load_genetic_code_overrides(gcode_dir: Path, genetic_code_id: int):
    """
    Dynamically loads codon-to-amino-acid overrides from tRNAscan-SE .gcode files.
    Skips comment lines starting with '#' and builds a custom translation map.
    """
    custom_code = STANDARD_GENETIC_CODE.copy()
    gcode_file = gcode_dir / f"{genetic_code_id}.gcode"
    
    if not gcode_file.exists():
        # Fallback safeguard for Candida albicans if the physical file is missing
        if genetic_code_id == 12:
            custom_code["CTG"] = "S"
        return custom_code
        
    with open(gcode_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            # Expected line format: <Codon> <3-letter AA> <1-letter AA> (e.g., TGA Trp W)
            parts = line.split()
            if len(parts) >= 3:
                codon = parts[0].upper()
                aa_one_letter = parts[2].upper()
                
                if len(codon) == 3 and len(aa_one_letter) == 1:
                    custom_code[codon] = aa_one_letter
                    
    return custom_code

def translate_cds(dna_seq, genetic_map):
    """Translates nucleotide sequence into amino acids using the provided mapping dictionary."""
    aa_seq = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        if len(codon) < 3:
            continue
        aa_seq.append(genetic_map.get(codon, "X"))
    return "".join(aa_seq)

def calculate_lcr_features(aa_seq, window_size=15, entropy_threshold=2.2):
    """Calculates LCR metrics and spatial localizations across 3 and 10 bins using Shannon Entropy."""
    lcr_regions = []
    L = len(aa_seq)
    if L < window_size:
        return 0, 0, [0,0,0], [0]*10, "", "", ""
    
    i = 0
    while i <= L - window_size:
        window = aa_seq[i:i+window_size]
        counts = Counter(window)
        entropy = -sum((c/window_size) * math.log2(c/window_size) for c in counts.values())
        
        if entropy < entropy_threshold:
            start = i
            while i <= L - window_size:
                w_next = aa_seq[i:i+window_size]
                c_next = Counter(w_next)
                e_next = -sum((c/window_size) * math.log2(c/window_size) for c in c_next.values())
                if e_next >= entropy_threshold:
                    break
                i += 1
            end = i + window_size - 1
            lcr_regions.append((start, end, aa_seq[start:end+1]))
        i += 1

    total_lcr_count = len(lcr_regions)
    total_lcr_length = sum(end - start + 1 for start, end, _ in lcr_regions)
    
    bins_3 = [0, 0, 0]
    bins_10 = [0] * 10
    
    n_seqs, m_seqs, c_seqs = [], [], []
    
    for start, end, seq in lcr_regions:
        rel_pos = ((start + end) / 2) / L
        
        if rel_pos <= 0.25:
            bins_3[0] += (end - start + 1)
            n_seqs.append(seq)
        elif rel_pos <= 0.75:
            bins_3[1] += (end - start + 1)
            m_seqs.append(seq)
        else:
            bins_3[2] += (end - start + 1)
            c_seqs.append(seq)
            
        bin_idx = min(int(rel_pos * 10), 9)
        bins_10[bin_idx] += (end - start + 1)
        
    return (
        total_lcr_count,
        total_lcr_length,
        bins_3,
        bins_10,
        ",".join(n_seqs) if n_seqs else "",
        ",".join(m_seqs) if m_seqs else "",
        ",".join(c_seqs) if c_seqs else ""
    )

def parse_cds_fasta(fasta_path, genetic_map):
    """Parses multi-FASTA CDS to collect sequence lengths and calculate LCR profiles using custom code mapping."""
    protein_data = {}
    current_id = None
    current_seq = []
    
    if not fasta_path.exists():
        return protein_data

    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    dna = "".join(current_seq)
                    aa = translate_cds(dna, genetic_map)
                    protein_data[current_id] = (len(aa), calculate_lcr_features(aa))
                current_id = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            dna = "".join(current_seq)
            aa = translate_cds(dna, genetic_map)
            protein_data[current_id] = (len(aa), calculate_lcr_features(aa))
            
    return protein_data

# --------------------
# --- Main function --
# --------------------
def main():
    args = arguments()
    
    dataset_df = pd.read_csv(args.dataset_sheet, sep="\t")
    active_samples = dataset_df[dataset_df['include'] == True].copy()
    
    compiled_rows = []
    
    for _, row in active_samples.iterrows():
        sample_id = row['sample']
        accession_id = row['accession']
        gcode = int(row['genetic_code'])
        
        # Dynamically load custom codon modifications from the .gcode file
        genetic_map = load_genetic_code_overrides(args.gcode_dir, gcode)
        
        summary_path = args.per_genome_dir / sample_id / "codon_metrics" / f"{sample_id}_summary.tsv"
        
        cds_pattern = str(args.cds_dir / row['cds'])
        from glob import glob
        matched_cds = glob(cds_pattern)
        
        if not summary_path.exists() or not matched_cds:
            print(f"Warning: Skipping {sample_id}. Missing summary metrics or CDS file context.")
            continue
            
        # Pass the dynamic genetic_map for translations and LCR calculations
        protein_structural_map = parse_cds_fasta(Path(matched_cds[0]), genetic_map)
        
        summary_df = pd.read_csv(summary_path, sep="\t")
        
        for _, gen_row in summary_df.iterrows():
            gene_id = gen_row['seqid'] if 'seqid' in summary_df.columns else gen_row.iloc[0]
            
            if gene_id in protein_structural_map:
                prot_len, lcr_meta = protein_structural_map[gene_id]
                lcr_count, lcr_len, b3, b10, n_seq, m_seq, c_seq = lcr_meta
            else:
                prot_len, lcr_count, lcr_len = 0, 0, 0
                b3, b10 = [0]*3, [0]*10
                n_seq, m_seq, c_seq = "", "", ""
                
            compiled_rows.append({
                "Assembly ID from NCBI": accession_id,
                "sample": sample_id,
                "Protein ID (NCBI accession)": gene_id,
                "Protein length": prot_len,
                "Presence of protein domains; Boolean": False,
                "symbolic localization of protein domains; 10 bins scaled to sum up to total protein length": ",".join(["0"]*10),
                "number of transmembrane elements predicted with TMHMM": 0,
                "total length of transmembrane elements": 0,
                "Symbolic localization of transmembrane elements; 10 bins scaled to sum up to total protein length": ",".join(["0"]*10),
                "Presence of signal peptide; Boolean": False,
                "Total number of LCR": lcr_count,
                "Total length of LCR": lcr_len,
                "Symbolic localization of LCR; 3 bins: N-termini (0-0.25 of protein length), middle (0.25-0.75 of protein length), and C-term (0.75-1 protein length)": ",".join(map(str, b3)),
                "Symbolic localization of LCR; 10 bins scaled to sum up to total protein length": ",".join(map(str, b10)),
                "LCR sequences in the N-terminal part of protein, separated by a comma": n_seq,
                "LCR sequences in the middle part of protein, separated by a comma": m_seq,
                "LCR sequences in the  C-terminal part of the protein, separated by a comma": c_seq,
                "Pfam domains overlapping with LCRs (>80% of LCR length)": "",
                "Pfam domains in protein (ordered by domain start)": "",
                "GO terms based on Pfam domains obtained by mapping on pfam2go, separated with the pipe symbol '|'": "",
                "kingdom": row.get('kingdom', 'Fungi'),
                "phylum": row.get('phylum', 'Unknown'),
                "lifestyle": row.get('lifestyle', 'Unknown'),
                "tAI": gen_row.get('tAI', min(gen_row.filter(like='tai').values, default=0.0)),
                "ENC": gen_row.get('ENC', min(gen_row.filter(like='enc').values, default=61.0)),
                "GC": gen_row.get('GC', min(gen_row.filter(like='gc').values, default=0.0)),
                "GC3s": gen_row.get('GC3s', min(gen_row.filter(like='gc3s').values, default=0.0))
            })
            
    out_df = pd.DataFrame(compiled_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"Success! Master table compiled with {len(out_df)} gene records saved to: {args.output}")

if __name__ == "__main__":
    main()