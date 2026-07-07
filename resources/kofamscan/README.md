KOfam: Profile HMM database of KEGG Orthology
====================================================

Laboratory of Chemical Life Science, Bioinformatics Center,
Institute for Chemical Research, Kyoto University
Uji, Kyoto 611-0011, Japan
Web site: https://www.genome.jp/tools/kofamkoala/
Feedback: https://www.genome.jp/feedback/
==================================================================

KOfam is a profile hidden Markov model (HMM) database of KEGG Orthology (KO).
Profiles are built from sequences in KO database using CD-HIT, MAFFT and
HMMER. Each profile has its own HMMER score threshold, with which the KO is
assigned to a sequence.
Annotation with KOfam can be done using a command line tool named KofamScan
or a web tool named KofamKOALA.

KofamScan:  ftp://ftp.genome.jp/pub/tools/kofamscan/
KofamKOALA: https://www.genome.jp/tools/kofamkoala/

This directory contains the following files.
profiles.tar.gz: profile HMMs of all KO groups.
ko_list.gz: Tab separated file containing the following information:
    knum         ... K number
    threshold    ... score threshold
    score_type   ... score type used for the KO (full or domain)
    profile_type ... Poorly aligned sequences are removed (trim) or not (all)
                     when building the profile
    F-measure    ... F-measure when calculating threshold
    nseq         ... number of sequences
    nseq_used    ... number of sequences used to build the profile
    alen         ... alignment length
    mlen         ... length of consensus positions
    eff_nseq     ... effective number of sequences
    re/pos       ... relative entropy per position
    definition   ... KO definition

==================================================================
Last update: 2019/4/9