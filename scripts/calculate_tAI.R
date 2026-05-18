# To calculate tAI:

# Usage
# get_tai(cf, trna_w, w_format = "cubar")
# Arguments
# cf A matrix of codon frequencies as calculated by count_codons(). Note: Start
# codons should be removed from sequences before analysis to avoid bias from
# universal start codon usage.
# trna_w A table of tRNA weights for each codon, generated using est_trna_weight().
# These weights reflect relative tRNA availability.
# w_format Character string specifying the format of tRNA weights: "cubar" (default, weights
# from cubar package) or "tAI" (weights from the tAI package format)






# instruction:
# https://cran.r-project.org/web/packages/cubar/cubar.pdf