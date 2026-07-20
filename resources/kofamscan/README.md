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



# KOfam resources for the ribosomal CAI reference

This directory contains the local KOfam metadata and the restricted reference
definition used to identify cytosolic eukaryotic ribosomal proteins. The
resulting CDS set is the sample-specific reference from which CAI codon weights
are estimated.

The workflow does not run an unrestricted whole-proteome KO annotation for
downstream functional analysis. It restricts KofamScan to the large- and
small-subunit cytosolic eukaryotic ribosomal KOs extracted from KEGG BRITE
`ko03011`.

## Required layout

```text
resources/kofamscan/
|-- ko_list
|-- profiles/
|   `-- Kxxxxx.hmm
|-- kofam.config.yaml
`-- ribosome_reference/
    |-- ko03011.json
    |-- SOURCE.yaml
    |-- cytosolic_eukaryotic_ribosome_kos.tsv
    `-- cytosolic_eukaryotic_ribosome.hal
```

`ko_list` and the BRITE reference files are present in this repository.
`profiles/` is not bundled because the HMM database is large; it must be
populated from an official KOfam release before the Snakemake DAG can be built.
The configured path is `kofamscan.database.profiles_dir`.

Use the uncompressed `ko_list` filename expected by `config/config.yaml` and
unpack the profile archive so that individual `Kxxxxx.hmm` files are directly
available under `profiles/`. The `ko_list` and HMM profiles should come from the
same KOfam release.

## What each file does

| Path | Role |
| --- | --- |
| `ko_list` | Provides one KO's score threshold, score type, profile type and definition. Thresholds determine whether KofamScan assignments are significant. |
| `profiles/Kxxxxx.hmm` | Profile HMM searched by `exec_annotation`. Only profiles named in the restricted HAL file are used. |
| `ribosome_reference/ko03011.json` | Local KEGG BRITE snapshot from which the eukaryotic cytosolic large- and small-subunit KO branches are extracted. |
| `ribosome_reference/SOURCE.yaml` | Records the BRITE source, retrieval date, included branches and excluded bacterial, archaeal and organellar branches. |
| `cytosolic_eukaryotic_ribosome_kos.tsv` | Auditable table of selected KOs, subunits, thresholds and profile metadata. |
| `cytosolic_eukaryotic_ribosome.hal` | Restricted list of HMM paths consumed by KofamScan. It is generated from the BRITE snapshot, `ko_list` and available profiles. |
| `kofam.config.yaml` | Standalone Kofam path template; the Snakemake rule passes its own explicit profile, KO-list, CPU and output arguments. |

## How the workflow uses these resources

1. `prepare_kofam_ribosomal_reference` parses `ko03011.json`, keeps only
   cytosolic eukaryotic ribosomal protein branches, joins KOfam thresholds and
   verifies that every selected HMM exists.
2. `run_kofamscan_ribosome` searches each sample proteome against the restricted
   HAL profile list with `exec_annotation`.
3. `parse_kofamscan_ribosome` retains threshold-passing hits and maps protein
   IDs to gene and CDS IDs through `gene_protein_map.tsv`.
4. `select_ribosomal_reference_cds` ranks significant hits and, by default,
   keeps the best CDS per KO. The selected nucleotide sequences are passed to
   `calculate_tAI.R` as the CAI reference.

This design avoids defining ribosomal genes through free-text FASTA headers and
prevents mitochondrial, bacterial or archaeal ribosomal profiles from entering
the default fungal CAI reference.

## Validation and update policy

The workflow stops if the profile directory is absent, a selected KO lacks a
profile, required thresholds are missing, identifiers cannot be mapped, or
fewer than five valid ribosomal reference CDS remain. Falling below
`ribosomal_reference.min_reference_genes` produces metric QC status `WARN`; it
stops execution only when `trna_qc.mode` is `error`.

When updating KOfam, replace `ko_list` and `profiles/` together. When updating
the BRITE snapshot, also update `SOURCE.yaml` and regenerate the restricted KO
table and HAL file through `prepare_kofam_ribosomal_reference`. Changes to any
of these resources can change the CAI reference and therefore require CAI and
downstream analyses to be recomputed.

KOfam and KofamScan are distributed by the KEGG/GenomeNet Bioinformatics
Center: <https://www.genome.jp/tools/kofamkoala/>.
