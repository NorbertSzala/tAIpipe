Scripts

    - count_tAI.R - script from CUBAR R Package counting tAI - main script.
    
    - count_CAI.R - script from CUBAR R package counting CAI and other usefull codon optimization statistics

    - codonM	script in Pearl what count a codons in sequence (CDS seq needed)

    - count_codons.py -  script used to make count_codons.txt. It calculates how many codons were in sequence, translate it into DNA and sort it in TCAG order.

    - codonW and tRNAscanse- install manually

    - tRNAscanse_pipeline_tAI_updated.sh - instructions how to run tRNAscan-SE.
    - other_analysis_pipeline_tAI.sh - aggregates results from codonM count_codons.py and codonW with tAI.R script. Count tAI is newer version




Files:

    - codons_NC_bulk	codonW output. counted CAI frequency of codons in each sequence. File with _total ic summary to whole CDS
    
    - codons_NC_output  codonW ouptut. counted some parameters like: NC, GC3s, GC, L_aa for each sequence. NC means how many different codons were used to code protein in that sequence. GC3s - number of synonymous codons with G or C in 3rd position. Laa- total number of synonymous and non-synomous codons


In usefull directory, you can find some usefull scripts that will help create graphical analyze to dashboard