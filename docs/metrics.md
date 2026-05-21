# Metrics

`tAIpipe` calculates codon usage and translational adaptation metrics from CDS and tRNA gene data.

## tRNA Adaptation Index

The **tRNA Adaptation Index (tAI)** measures how well codons in a coding sequence match the available tRNA pool in a cell. It uses tRNA gene copy numbers as an approximation of tRNA availability and accounts for codon–anticodon pairing, including wobble interactions.

Higher tAI values indicate that a gene uses codons recognized by more abundant tRNAs. This may support faster or more efficient translation. Lower local tAI may also be biologically meaningful, because slower translation can contribute to protein folding or regulatory pauses.

## CAI

The **Codon Adaptation Index (CAI)** measures how similar a gene's codon usage is to a reference set of highly expressed genes. CAI values range from 0 to 1. Higher values indicate stronger similarity to preferred codon usage.

## ENC

The **Effective Number of Codons (ENC)** describes the overall strength of codon usage bias. Lower ENC values indicate stronger codon bias, while higher values indicate more uniform use of synonymous codons.

## RSCU

**Relative Synonymous Codon Usage (RSCU)** compares observed codon usage with expected codon usage under equal use of synonymous codons. RSCU > 1 indicates overrepresentation, while RSCU < 1 indicates underrepresentation.

## GC and GC3s

`GC` measures overall GC content. `GC3s` measures GC content at synonymous third codon positions and is useful for describing codon-level nucleotide bias.