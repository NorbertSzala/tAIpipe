# tAIpipe
Bioinformatics pipeline for tRNA adaptation index (tAI) analysis across CDS and genomic datasets, with an interactive dashboard.

## Goal

The main goal of this project is to create interactive dashboard integrating data from tRNA Adaptation Index analysis. User can choose features what are plotted, eg. tAI vs genome size. This pipeline does not limit to fungal species.

## Requirements
    - python
    - R
    - tRNAscan-SE
    - HMMscan


## Pipeline 

The first step is to define what organisms should be analyzed (by NCBI accession number or by species name. In second case there will be downloaded reference genome).

Then pipeline creates './data/metadata/samples.tsv' files with desc