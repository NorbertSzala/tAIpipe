# tAIpipe Analytical Report: Global and Local Translation Landscape

------------------------------------------------------------------------

## Quality Control & Dataset Auditing

### Cohort Ingestion Summary

#### Goal

Verify data ingestion and experiment layout integrity before
computational evaluation. Ensure taxonomic and ecological categories are
represented.

#### Methodology

Compile a cohort overview table summarizing total counts of independent
operational taxonomic units (samples) across Domains, Kingdoms, Phyla,
and Lifestyles.

#### Critical Observations

Look for extreme class imbalances or sample omissions. If a specific
lifestyle or phylum contains too few records compared to the rest of the
cohort, a text warning must alert the user to potential statistical
inflation during downstream grouping.

| domain  | kingdom | phylum        | lifestyle              | Total_Samples |
|:--------|:--------|:--------------|:-----------------------|--------------:|
| eukarya | Fungi   | Ascomycota    | nectar/tap_saprotroph  |             1 |
| eukarya | Fungi   | Ascomycota    | nectar_tap_saprotroph  |             1 |
| eukarya | Fungi   | Ascomycota    | plant_pathogen         |             1 |
| eukarya | Fungi   | Ascomycota    | unspecified_saprotroph |             1 |
| eukarya | Fungi   | Basidiomycota | wood_saprotroph        |             1 |

Summary of Successfully Ingested Cohort Groups

<strong>⚠️ WARNING:</strong> The following sample cohorts have fewer
than 3 active biological representatives. High risk of statistical
inflation during macro comparative groupings:

| domain  | kingdom | phylum        | lifestyle              | Total_Samples |
|:--------|:--------|:--------------|:-----------------------|--------------:|
| eukarya | Fungi   | Ascomycota    | nectar/tap_saprotroph  |             1 |
| eukarya | Fungi   | Ascomycota    | nectar_tap_saprotroph  |             1 |
| eukarya | Fungi   | Ascomycota    | plant_pathogen         |             1 |
| eukarya | Fungi   | Ascomycota    | unspecified_saprotroph |             1 |
| eukarya | Fungi   | Basidiomycota | wood_saprotroph        |             1 |

### Organism-Level Assembly & Prediction QC

#### Goal

Identify analytical outliers, structural sequence assembly fragmentation
artifacts, and potential downstream tRNA prediction errors.

#### Methodology

Generate a dual panel of sample-wide scatter plots: 1. Genome Assembly
Size (Mbp) vs. Total Predicted tRNA Gene Count. 2. Genome-wide Average
GC Content vs. Genome-wide Average Synonymous Third-Position GC Content
(GC3s).

#### Critical Observations

Legitimate biological trends exhibit linear scaling between genome size
and tRNA pools. Mark severe technical outliers in red with text labels.
For instance, a bloated genome sequence with an artificially suppressed
tRNA count indicates low assembly coverage or raw file contamination.

#### Notes

Samples identified as suspicious through the entire pipeline auditing
process (across relevant data splits) must be explicitly flagged and
listed here.

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/assembly-prediction-qc-plots-1.png)

------------------------------------------------------------------------

## Evolutionary & Cohort-Scale Analysis (Across All Species)

### Translational Strategy Mapped by Lifestyle and Taxonomy

#### Goal

Evaluate how broad macro-evolutionary descent (Phylum) and ecological
pressures (Lifestyle) dictate global codon usage selection across the
entire cohort landscape.

#### Methodology

Compile an integrated multi-panel visualization mapping genome-wide mean
tAI against mean Effective Number of Codons (ENC/Nc). Points represent
distinct species, color-coded by lifestyle categories and shaped by
taxonomic phyla, paired with marginal distribution boxplots.

#### Critical Observations

Determine if specific niche adaptations (e.g., high-impact plant
pathogens vs. free-living saprotrophs) cluster into separated spaces of
high translational efficiency (high tAI) paired with narrowed codon
usage (low ENC).

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/translational-strategy-macro-1.png)

### Selection Pressure vs Mutational Drift (The Wright’s Curve)

#### Goal

Disentangle whether observed codon usage bias across species is driven
by natural selection for translational speed/accuracy or by passive
genomic background mutational drift.

#### Methodology

Map genome-wide average ENC against average GC3s. Overlay the
theoretical mathematical Wright’s curve (*H*<sub>0</sub>: purely neutral
mutational drift) defined by:
$$ENC\_{\text{theoretical}} = 2 + s + \frac{29}{s^2 + (1-s)^2}$$
where *s* = *G**C*3*s*.

#### Critical Observations

Species or functional cohorts pooling significantly below the
theoretical line undergo active directional selection toward optimal
codons. Species hugging the curve evolve primarily under neutral
mutational drift.

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/wrights-curve-analysis-1.png)

------------------------------------------------------------------------

## Type-Level Functional Protein & Gene Adaptation (The Extremes)

### Structural Features of Protein Extremes Across Cohorts

#### Goal

Determine global amino acid and protein structural sorting rules
relative to translational efficiency across the entire experimental
framework.

#### Methodology

Aggregate gene sequences from all genomes into a pooled dataset. Isolate
top 10% (optimized) and bottom 10% (slowed/regulatory) global tAI
thresholds. Construct dodge bar charts displaying the fraction of genes
containing signal peptides or transmembrane domains (TM), grouped by
Phylum or Lifestyle.

#### Spatial Profiling

For the Top 10% and Bottom 10% tAI groups, generate density profile
charts across 10 unified length bins tracking the positional layout of:
\* Predicted protein domain boundaries (Pfam) \* Transmembrane segments
(TMHMM) \* Low Complexity Regions (LCR) This directly checks if highly
optimized proteins restrict sequence complexity or transmembrane blocks
to specific structural domains (e.g., N- or C-termini) compared to
slower translational targets.

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/cross-species-structural-extremes-1.png)

### Global Pfam & GO Terms Enrichment of Translation Extremes

#### Goal

Uncover universal metabolic or evolutionary pathways that are strictly
conserved for rapid translation (Top) or translational pausing (Bottom)
across taxonomic boundaries.

#### Methodology

Run hyper-geometric distribution tests for overrepresentation of Pfam
domain keys and Gene Ontology (GO) terms (derived via `pfam2go`
mappings) within the combined Top/Bottom 10% tAI gene pools.

#### Critical Observations

Results are presented as profiles of top enriched functions broken down
by Phylum and Lifestyle to highlight potential instances of evolutionary
convergence among ecologically matching pathogens.

| tAI_Group | phylum | GO_raw | Count_Selected | Universe_N | Sample_Size | p_value | adj_p_value |
|:--------|:------|:------|------------:|---------:|:----------|-------:|----------:|

Top 20 Broad Functional GO Terms Enriched in Translation Rate Extremes

------------------------------------------------------------------------

## Individual Organism Deep-Dive (Per-Genome Diagnostics)

*Note: This diagnostic loop evaluates target individual assemblies
dynamically or cycles sequentially through distinct species
configurations.*

### Ztritici

#### Intra-Organismal tAI Distribution and Core Outliers

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-organism-deep-dive-loop-1.png)
\##### Top 10 High-Efficiency Expression Target Genes (Elite Pool)

| Protein ID (NCBI accession) | tAI | ENC | GO_Terms |
|:--------------------------------------------------------|---:|--:|:--------|
| lcl\|CM001196.1_cds_EGP91094.1_572 \[locus_tag=MYCGRDRAFT_78335\] \[db_xref=JGIDB:Mycgr3_78335\] \[protein=hypothetical protein\] \[protein_id=EGP91094.1\] \[location=1906441..1906695\] \[gbkey=CDS\] | 0.8695396 | 50.04355 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91548.1_1446 \[locus_tag=MYCGRDRAFT_102767\] \[db_xref=JGIDB:Mycgr3_102767\] \[protein=hypothetical protein\] \[protein_id=EGP91548.1\] \[location=4368971..4369159\] \[gbkey=CDS\] | 0.8681674 | 47.83466 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91714.1_1779 \[locus_tag=MYCGRDRAFT_90143\] \[db_xref=JGIDB:Mycgr3_90143\] \[protein=hypothetical protein\] \[protein_id=EGP91714.1\] \[location=5396657..5396980\] \[gbkey=CDS\] | 0.8676202 | 49.40411 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91119.1_618 \[locus_tag=MYCGRDRAFT_102256\] \[db_xref=JGIDB:Mycgr3_102256\] \[protein=hypothetical protein\] \[protein_id=EGP91119.1\] \[location=2030498..2030665\] \[gbkey=CDS\] | 0.8673825 | 45.61720 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91273.1_913 \[locus_tag=MYCGRDRAFT_83763\] \[db_xref=JGIDB:Mycgr3_83763\] \[protein=hypothetical protein\] \[protein_id=EGP91273.1\] \[location=2856073..2856261\] \[gbkey=CDS\] | 0.8659387 | 52.94876 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP92619.1_356 \[locus_tag=MYCGRDRAFT_65284\] \[db_xref=JGIDB:Mycgr3_65284\] \[protein=hypothetical protein\] \[protein_id=EGP92619.1\] \[location=complement(join(1324674..1324948,1325030..1325366))\] \[gbkey=CDS\] | 0.8648925 | 53.19686 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91636.1_1626 \[locus_tag=MYCGRDRAFT_102855\] \[db_xref=JGIDB:Mycgr3_102855\] \[protein=hypothetical protein\] \[protein_id=EGP91636.1\] \[location=4892793..4893053\] \[gbkey=CDS\] | 0.8643951 | 52.84462 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP92169.1_1264 \[locus_tag=MYCGRDRAFT_66663\] \[db_xref=InterPro:IPR007109,JGIDB:Mycgr3_66663\] \[protein=hypothetical protein\] \[frame=2\] \[partial=5’\] \[protein_id=EGP92169.1\] \[location=complement(3861310..\>3862198)\] \[gbkey=CDS\] | 0.8640481 | 51.02636 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91550.1_1448 \[locus_tag=MYCGRDRAFT_107522\] \[db_xref=JGIDB:Mycgr3_107522\] \[protein=hypothetical protein\] \[protein_id=EGP91550.1\] \[location=join(4374013..4374077,4374181..4374294,4374778..4375100,4375157..4375410,4375460..4375798)\] \[gbkey=CDS\] | 0.8637892 | 49.20338 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91429.1_1229 \[locus_tag=MYCGRDRAFT_102643\] \[db_xref=JGIDB:Mycgr3_102643\] \[protein=hypothetical protein\] \[protein_id=EGP91429.1\] \[location=3739650..3739931\] \[gbkey=CDS\] | 0.8636675 | 52.06929 | Unannotated / No Pfam Domain Found |

##### Bottom 10 Regulatory/Slowing Translation Genes (Structural Control)

| Protein ID (NCBI accession) | tAI | ENC | GO_Terms |
|:--------------------------------------------------------|---:|--:|:--------|
| lcl\|CM001196.1_cds_EGP91819.1_1984 \[locus_tag=MYCGRDRAFT_79161\] \[db_xref=JGIDB:Mycgr3_79161\] \[protein=hypothetical protein\] \[protein_id=EGP91819.1\] \[location=complement(join(6006196..6006264,6006336..6006407,6006494..6006559))\] \[gbkey=CDS\] | 0.8234634 | 49.13446 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP92200.1_1195 \[locus_tag=MYCGRDRAFT_89647\] \[db_xref=JGIDB:Mycgr3_89647\] \[protein=hypothetical protein\] \[protein_id=EGP92200.1\] \[location=complement(join(3644229..3644323,3644388..3644474,3644596..3644737,3644793..3644798))\] \[gbkey=CDS\] | 0.8268069 | 51.46783 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91658.1_1661 \[locus_tag=MYCGRDRAFT_84191\] \[db_xref=JGIDB:Mycgr3_84191\] \[protein=hypothetical protein\] \[protein_id=EGP91658.1\] \[location=5026169..5026330\] \[gbkey=CDS\] | 0.8294889 | 53.10237 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91577.1_1501 \[locus_tag=MYCGRDRAFT_102797\] \[db_xref=InterPro:IPR006461,JGIDB:Mycgr3_102797\] \[protein=hypothetical protein\] \[protein_id=EGP91577.1\] \[location=join(4545681..4545764,4545818..4546252)\] \[gbkey=CDS\] | 0.8301620 | 47.33881 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91953.1_1716 \[locus_tag=MYCGRDRAFT_90092\] \[db_xref=JGIDB:Mycgr3_90092\] \[protein=hypothetical protein\] \[protein_id=EGP91953.1\] \[location=complement(5235360..5235686)\] \[gbkey=CDS\] | 0.8314079 | 48.61781 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91585.1_1519 \[locus_tag=MYCGRDRAFT_53984\] \[db_xref=InterPro:IPR003213,JGIDB:Mycgr3_53984\] \[protein=hypothetical protein\] \[protein_id=EGP91585.1\] \[location=join(4590613..4590661,4590821..4590936,4590993..4591091)\] \[gbkey=CDS\] | 0.8318188 | 45.80838 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP91727.1_1804 \[locus_tag=MYCGRDRAFT_102916\] \[db_xref=JGIDB:Mycgr3_102916\] \[protein=hypothetical protein\] \[protein_id=EGP91727.1\] \[location=5459538..5459837\] \[gbkey=CDS\] | 0.8332747 | 51.11220 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP90871.1_129 \[locus_tag=MYCGRDRAFT_106780\] \[db_xref=JGIDB:Mycgr3_106780\] \[protein=hypothetical protein\] \[protein_id=EGP90871.1\] \[location=577551..578981\] \[gbkey=CDS\] | 0.8345274 | 40.77112 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP90890.1_170 \[locus_tag=MYCGRDRAFT_102033\] \[db_xref=JGIDB:Mycgr3_102033\] \[protein=hypothetical protein\] \[protein_id=EGP90890.1\] \[location=723354..723716\] \[gbkey=CDS\] | 0.8359453 | 51.99785 | Unannotated / No Pfam Domain Found |
| lcl\|CM001196.1_cds_EGP90865.1_118 \[locus_tag=MYCGRDRAFT_88698\] \[db_xref=JGIDB:Mycgr3_88698\] \[protein=hypothetical protein\] \[protein_id=EGP90865.1\] \[location=join(536720..536772,536829..537015,537089..537157,537226..537396)\] \[gbkey=CDS\] | 0.8364195 | 46.15146 | Unannotated / No Pfam Domain Found |

#### Local Protein Feature Correlation

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-organism-deep-dive-loop-2.png)
\#### Overlap of Structural Domains and Sequence Complexity

| Protein ID (NCBI accession) | Local_Cohort | Overlapping_Domains |
|:------------------------------------------------------------|:----|:-----|
| lcl\|CM001196.1_cds_EGP90808.1_5 \[locus_tag=MYCGRDRAFT_51803\] \[db_xref=InterPro:IPR000571,JGIDB:Mycgr3_51803\] \[protein=hypothetical protein\] \[protein_id=EGP90808.1\] \[location=120041..120670\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP92785.1_10 \[locus_tag=MYCGRDRAFT_88592\] \[db_xref=InterPro:IPR007867,JGIDB:Mycgr3_88592\] \[protein=hypothetical protein\] \[protein_id=EGP92785.1\] \[location=complement(join(132334..132494,132758..133296,133366..133440,133734..133763,133968..134113))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP90811.1_12 \[locus_tag=MYCGRDRAFT_101983\] \[db_xref=JGIDB:Mycgr3_101983\] \[protein=hypothetical protein\] \[protein_id=EGP90811.1\] \[location=136333..136869\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP92781.1_15 \[locus_tag=MYCGRDRAFT_88597\] \[db_xref=JGIDB:Mycgr3_88597\] \[protein=hypothetical protein\] \[protein_id=EGP92781.1\] \[location=complement(145496..146929)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP90818.1_31 \[locus_tag=MYCGRDRAFT_88613\] \[db_xref=JGIDB:Mycgr3_88613\] \[protein=hypothetical protein\] \[protein_id=EGP90818.1\] \[location=join(206531..206728,206987..207093,207754..207760,207914..208076,208658..208710,209351..209390,209681..209717,209998..210041,210291..210351,210603..210654,210916..211064,211149..211224)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP92771.1_32 \[locus_tag=MYCGRDRAFT_88614\] \[db_xref=JGIDB:Mycgr3_88614\] \[protein=hypothetical protein\] \[protein_id=EGP92771.1\] \[location=complement(211916..212437)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP90819.1_34 \[locus_tag=MYCGRDRAFT_64883\] \[db_xref=JGIDB:Mycgr3_64883\] \[protein=hypothetical protein\] \[protein_id=EGP90819.1\] \[location=join(217709..218075,218241..218599)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP90822.1_39 \[locus_tag=MYCGRDRAFT_88619\] \[db_xref=JGIDB:Mycgr3_88619\] \[protein=hypothetical protein\] \[protein_id=EGP90822.1\] \[location=join(234190..234686,234995..235022)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP90823.1_41 \[locus_tag=MYCGRDRAFT_106743\] \[db_xref=JGIDB:Mycgr3_106743\] \[protein=hypothetical protein\] \[protein_id=EGP90823.1\] \[location=237862..238635\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP92766.1_42 \[locus_tag=MYCGRDRAFT_88622\] \[db_xref=JGIDB:Mycgr3_88622\] \[protein=hypothetical protein\] \[protein_id=EGP92766.1\] \[location=complement(239195..240148)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP90824.1_44 \[locus_tag=MYCGRDRAFT_88624\] \[db_xref=JGIDB:Mycgr3_88624\] \[protein=hypothetical protein\] \[protein_id=EGP90824.1\] \[location=join(246303..246499,246564..246675)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP92763.1_49 \[locus_tag=MYCGRDRAFT_88629\] \[db_xref=JGIDB:Mycgr3_88629\] \[protein=hypothetical protein\] \[protein_id=EGP92763.1\] \[location=complement(265655..266656)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP92762.1_50 \[locus_tag=MYCGRDRAFT_88630\] \[db_xref=JGIDB:Mycgr3_88630\] \[protein=hypothetical protein\] \[protein_id=EGP92762.1\] \[location=complement(269253..270098)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP92757.1_64 \[locus_tag=MYCGRDRAFT_88643\] \[db_xref=JGIDB:Mycgr3_88643\] \[protein=hypothetical protein\] \[protein_id=EGP92757.1\] \[location=complement(join(313055..313083,313373..313620,313881..313957))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP92752.1_70 \[locus_tag=MYCGRDRAFT_88648\] \[db_xref=JGIDB:Mycgr3_88648\] \[protein=hypothetical protein\] \[protein_id=EGP92752.1\] \[location=complement(join(345834..346355,347133..347144))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP90840.1_73 \[locus_tag=MYCGRDRAFT_106752\] \[db_xref=InterPro:IPR002110,InterPro:IPR007087,JGIDB:Mycgr3_106752\] \[protein=hypothetical protein\] \[protein_id=EGP90840.1\] \[location=353858..356317\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP92750.1_75 \[locus_tag=MYCGRDRAFT_32102\] \[db_xref=InterPro:IPR002734,JGIDB:Mycgr3_32102\] \[protein=hypothetical protein\] \[protein_id=EGP92750.1\] \[location=complement(358191..358973)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP92749.1_78 \[locus_tag=MYCGRDRAFT_64923\] \[db_xref=InterPro:IPR002198,JGIDB:Mycgr3_64923\] \[protein=hypothetical protein\] \[protein_id=EGP92749.1\] \[location=complement(join(375594..375842,375909..376460))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP90844.1_81 \[locus_tag=MYCGRDRAFT_101995\] \[db_xref=InterPro:IPR001841,JGIDB:Mycgr3_101995\] \[protein=hypothetical protein\] \[protein_id=EGP90844.1\] \[location=383038..383349\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CM001196.1_cds_EGP90845.1_82 \[locus_tag=MYCGRDRAFT_35093\] \[db_xref=InterPro:IPR006811,JGIDB:Mycgr3_35093\] \[protein=hypothetical protein\] \[protein_id=EGP90845.1\] \[location=join(384317..384474,384528..384807,384871..385260)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |

### Pplacenta

#### Intra-Organismal tAI Distribution and Core Outliers

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-organism-deep-dive-loop-3.png)
\##### Top 10 High-Efficiency Expression Target Genes (Elite Pool)

| Protein ID (NCBI accession) | tAI | ENC | GO_Terms |
|:----------------------------------------------------------|--:|--:|:-------|
| lcl\|EQ966233.1_cds_EED86010.1_149 \[locus_tag=POSPLDRAFT_91540\] \[protein=predicted protein\] \[protein_id=EED86010.1\] \[location=join(1436785..1437049,1437106..1437198,1437260..1437309)\] \[gbkey=CDS\] | 0.6457431 | 40.53922 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED85986.1_102 \[locus_tag=POSPLDRAFT_87690\] \[protein=predicted protein\] \[protein_id=EED85986.1\] \[location=join(969800..970252,970286..970336,971028..971586,971638..971783)\] \[gbkey=CDS\] | 0.6433874 | 43.52229 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED86024.1_24 \[locus_tag=POSPLDRAFT_91307\] \[protein=predicted protein\] \[protein_id=EED86024.1\] \[location=complement(227589..228212)\] \[gbkey=CDS\] | 0.6413327 | 45.68170 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED86066.1_104 \[locus_tag=POSPLDRAFT_91446\] \[protein=predicted protein\] \[protein_id=EED86066.1\] \[location=complement(join(982246..982489,982524..982627,982651..982974))\] \[gbkey=CDS\] | 0.6381110 | 49.61873 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED85989.1_109 \[locus_tag=POSPLDRAFT_91453\] \[db_xref=InterPro:IPR006016\] \[protein=predicted protein\] \[protein_id=EED85989.1\] \[location=join(1022774..1022815,1022880..1022918,1023229..1024118,1024172..1024206,1024268..1024326,1024378..1024498,1024562..1024628,1024679..1024736,1024794..1025057)\] \[gbkey=CDS\] | 0.6375067 | 42.14289 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED86074.1_119 \[locus_tag=POSPLDRAFT_108796\] \[db_xref=InterPro:IPR005197\] \[protein=hypothetical protein\] \[protein_id=EED86074.1\] \[location=complement(join(1119402..1119659,1119725..1119879,1119942..1120368,1120429..1120768,1121035..1121288))\] \[gbkey=CDS\] | 0.6367094 | 45.03399 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED86068.1_108 \[locus_tag=POSPLDRAFT_91452\] \[protein=predicted protein\] \[protein_id=EED86068.1\] \[location=complement(join(1020560..1020667,1020723..1021848,1021973..1022057,1022104..1022119))\] \[gbkey=CDS\] | 0.6345121 | 44.75709 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED86076.1_122 \[locus_tag=POSPLDRAFT_91476\] \[protein=predicted protein\] \[protein_id=EED86076.1\] \[location=complement(join(1142950..1143251,1143333..1143552,1143657..1143704,1143878..1143954,1144012..1144205,1144302..1144470,1144523..1144565))\] \[gbkey=CDS\] | 0.6334117 | 53.99999 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED86030.1_32 \[locus_tag=POSPLDRAFT_91322\] \[protein=predicted protein\] \[protein_id=EED86030.1\] \[location=complement(join(281403..281599,281659..281780,281832..282352,282403..282603))\] \[gbkey=CDS\] | 0.6328900 | 50.23630 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED85936.1_1 \[locus_tag=POSPLDRAFT_91257\] \[protein=predicted protein\] \[protein_id=EED85936.1\] \[location=join(5935..6007,6095..6372,6430..6592,6626..7251)\] \[gbkey=CDS\] | 0.6297691 | 51.75511 | Unannotated / No Pfam Domain Found |

##### Bottom 10 Regulatory/Slowing Translation Genes (Structural Control)

| Protein ID (NCBI accession) | tAI | ENC | GO_Terms |
|:------------------------------------------------------------|--:|--:|:-----|
| lcl\|EQ966233.1_cds_EED85947.1_22 \[locus_tag=POSPLDRAFT_91303\] \[protein=predicted protein\] \[protein_id=EED85947.1\] \[location=join(170275..170298,170665..170759,170815..171083,171303..171436,171493..171664,171714..171949)\] \[gbkey=CDS\] | 0.5793809 | 49.38096 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED85960.1_46 \[locus_tag=POSPLDRAFT_91343\] \[db_xref=InterPro:IPR000626\] \[protein=predicted protein\] \[protein_id=EED85960.1\] \[location=join(380259..380285,380542..380624,380682..380779,380832..381167,381224..381514,381578..381612)\] \[gbkey=CDS\] | 0.5864210 | 55.73391 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED86043.1_62 \[locus_tag=POSPLDRAFT_91374\] \[protein=predicted protein\] \[protein_id=EED86043.1\] \[location=complement(join(558055..558091,558222..558323,558408..558432,558483..558552,558602..558637,558679..558760,558793..558893,558944..558986,559042..559057,559115..559312,559361..559419,559449..559474,559539..559566,559678..559739,559790..559826,559882..559895,559962..560040,560098..560168))\] \[gbkey=CDS\] | 0.5870494 | 56.18197 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED85967.1_61 \[locus_tag=POSPLDRAFT_91371\] \[protein=predicted protein\] \[protein_id=EED85967.1\] \[location=join(543168..543191,543490..543516,543558..543644,543922..544062,544117..544142,544199..544224,544275..544318,544587..544640,544777..544797,545047..545116,545169..545193,545253..545412)\] \[gbkey=CDS\] | 0.5875857 | 53.32891 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED86036.1_52 \[locus_tag=POSPLDRAFT_91353\] \[protein=predicted protein\] \[protein_id=EED86036.1\] \[location=complement(join(439151..439203,439266..439414,439480..439504,439554..439615,439752..439811,440009..440098,440129..440144,440357..440384,440451..440476,440687..440766,440863..440966,441027..441105,441155..441222))\] \[gbkey=CDS\] | 0.5881783 | 56.99566 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED85968.1_64 \[locus_tag=POSPLDRAFT_91377\] \[protein=predicted protein\] \[protein_id=EED85968.1\] \[location=join(593877..593947,594003..594055,594111..594124,594200..594201,594257..594351,594410..594487,594538..594700,594755..594887)\] \[gbkey=CDS\] | 0.5919621 | 46.23928 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED85990.1_111 \[locus_tag=POSPLDRAFT_43667\] \[protein=hypothetical protein\] \[protein_id=EED85990.1\] \[location=join(1040896..1041094,1041173..1041477,1041532..1042185)\] \[gbkey=CDS\] | 0.5924709 | 57.25534 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED85991.1_113 \[locus_tag=POSPLDRAFT_87698\] \[db_xref=InterPro:IPR002893\] \[protein=predicted protein\] \[protein_id=EED85991.1\] \[location=join(1052030..1052358,1052436..1052564,1052638..1052642,1052694..1053295)\] \[gbkey=CDS\] | 0.5928463 | 54.52791 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED86009.1_148 \[locus_tag=POSPLDRAFT_91538\] \[protein=predicted protein\] \[protein_id=EED86009.1\] \[location=join(1426802..1426848,1426891..1427098,1427141..1427443)\] \[gbkey=CDS\] | 0.5930575 | 54.97698 | Unannotated / No Pfam Domain Found |
| lcl\|EQ966233.1_cds_EED86087.1_147 \[locus_tag=POSPLDRAFT_91537\] \[protein=predicted protein\] \[protein_id=EED86087.1\] \[location=complement(join(1423094..1423459,1423516..1423652,1423714..1423718,1423774..1423877,1423931..1424153,1424236..1424253,1424288..1424361,1424430..1424498,1424671..1424880,1424942..1424948,1425001..1425119,1425174..1425224))\] \[gbkey=CDS\] | 0.5933788 | 54.16997 | Unannotated / No Pfam Domain Found |

#### Local Protein Feature Correlation

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-organism-deep-dive-loop-4.png)
\#### Overlap of Structural Domains and Sequence Complexity

| Protein ID (NCBI accession) | Local_Cohort | Overlapping_Domains |
|:-------------------------------------------------------------|:----|:----|
| lcl\|EQ966233.1_cds_EED85936.1_1 \[locus_tag=POSPLDRAFT_91257\] \[protein=predicted protein\] \[protein_id=EED85936.1\] \[location=join(5935..6007,6095..6372,6430..6592,6626..7251)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED86015.1_3 \[locus_tag=POSPLDRAFT_91262\] \[protein=predicted protein\] \[protein_id=EED86015.1\] \[location=complement(join(22030..22073,22132..22243,22288..22332,22552..22898,22939..22957,23018..23338))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85941.1_12 \[locus_tag=POSPLDRAFT_91274\] \[protein=predicted protein\] \[protein_id=EED85941.1\] \[location=join(72283..72616,72709..73322)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85947.1_22 \[locus_tag=POSPLDRAFT_91303\] \[protein=predicted protein\] \[protein_id=EED85947.1\] \[location=join(170275..170298,170665..170759,170815..171083,171303..171436,171493..171664,171714..171949)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED86024.1_24 \[locus_tag=POSPLDRAFT_91307\] \[protein=predicted protein\] \[protein_id=EED86024.1\] \[location=complement(227589..228212)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED86030.1_32 \[locus_tag=POSPLDRAFT_91322\] \[protein=predicted protein\] \[protein_id=EED86030.1\] \[location=complement(join(281403..281599,281659..281780,281832..282352,282403..282603))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85959.1_44 \[locus_tag=POSPLDRAFT_91341\] \[protein=predicted protein\] \[protein_id=EED85959.1\] \[location=join(372276..372304,372420..372496,372566..372621,372674..372745,373045..373146,373203..373291,373345..373810)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85960.1_46 \[locus_tag=POSPLDRAFT_91343\] \[db_xref=InterPro:IPR000626\] \[protein=predicted protein\] \[protein_id=EED85960.1\] \[location=join(380259..380285,380542..380624,380682..380779,380832..381167,381224..381514,381578..381612)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED86036.1_52 \[locus_tag=POSPLDRAFT_91353\] \[protein=predicted protein\] \[protein_id=EED86036.1\] \[location=complement(join(439151..439203,439266..439414,439480..439504,439554..439615,439752..439811,440009..440098,440129..440144,440357..440384,440451..440476,440687..440766,440863..440966,441027..441105,441155..441222))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85967.1_61 \[locus_tag=POSPLDRAFT_91371\] \[protein=predicted protein\] \[protein_id=EED85967.1\] \[location=join(543168..543191,543490..543516,543558..543644,543922..544062,544117..544142,544199..544224,544275..544318,544587..544640,544777..544797,545047..545116,545169..545193,545253..545412)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED86043.1_62 \[locus_tag=POSPLDRAFT_91374\] \[protein=predicted protein\] \[protein_id=EED86043.1\] \[location=complement(join(558055..558091,558222..558323,558408..558432,558483..558552,558602..558637,558679..558760,558793..558893,558944..558986,559042..559057,559115..559312,559361..559419,559449..559474,559539..559566,559678..559739,559790..559826,559882..559895,559962..560040,560098..560168))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85968.1_64 \[locus_tag=POSPLDRAFT_91377\] \[protein=predicted protein\] \[protein_id=EED85968.1\] \[location=join(593877..593947,594003..594055,594111..594124,594200..594201,594257..594351,594410..594487,594538..594700,594755..594887)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85986.1_102 \[locus_tag=POSPLDRAFT_87690\] \[protein=predicted protein\] \[protein_id=EED85986.1\] \[location=join(969800..970252,970286..970336,971028..971586,971638..971783)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED86065.1_103 \[locus_tag=POSPLDRAFT_91445\] \[db_xref=InterPro:IPR000038\] \[protein=predicted protein\] \[protein_id=EED86065.1\] \[location=complement(join(978875..978888,978957..979071,979137..980127,980180..980226,980380..980691,980743..981129))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED86066.1_104 \[locus_tag=POSPLDRAFT_91446\] \[protein=predicted protein\] \[protein_id=EED86066.1\] \[location=complement(join(982246..982489,982524..982627,982651..982974))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED86068.1_108 \[locus_tag=POSPLDRAFT_91452\] \[protein=predicted protein\] \[protein_id=EED86068.1\] \[location=complement(join(1020560..1020667,1020723..1021848,1021973..1022057,1022104..1022119))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85989.1_109 \[locus_tag=POSPLDRAFT_91453\] \[db_xref=InterPro:IPR006016\] \[protein=predicted protein\] \[protein_id=EED85989.1\] \[location=join(1022774..1022815,1022880..1022918,1023229..1024118,1024172..1024206,1024268..1024326,1024378..1024498,1024562..1024628,1024679..1024736,1024794..1025057)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85990.1_111 \[locus_tag=POSPLDRAFT_43667\] \[protein=hypothetical protein\] \[protein_id=EED85990.1\] \[location=join(1040896..1041094,1041173..1041477,1041532..1042185)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85991.1_113 \[locus_tag=POSPLDRAFT_87698\] \[db_xref=InterPro:IPR002893\] \[protein=predicted protein\] \[protein_id=EED85991.1\] \[location=join(1052030..1052358,1052436..1052564,1052638..1052642,1052694..1053295)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|EQ966233.1_cds_EED85992.1_114 \[locus_tag=POSPLDRAFT_91462\] \[db_xref=InterPro:IPR005162\] \[protein=predicted protein\] \[protein_id=EED85992.1\] \[location=join(1061298..1062123,1062202..1062278)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |

### Spombe

#### Intra-Organismal tAI Distribution and Core Outliers

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-organism-deep-dive-loop-5.png)
\##### Top 10 High-Efficiency Expression Target Genes (Elite Pool)

| Protein ID (NCBI accession) | tAI | ENC | GO_Terms |
|:------------------------------------------------------------|--:|--:|:-----|
| lcl\|CU329670.1_cds_CAA93291.2_1822 \[gene=cdc8\] \[locus_tag=SPOM_SPAC27F1.02C\] \[db_xref=EnsemblGenomes-Gn:SPAC27F1.02c,EnsemblGenomes-Tr:SPAC27F1.02c.1,GOA:Q02088,InterPro:IPR000533,PomBase:SPAC27F1.02c,PomBase:SPAC27F1.02c.1,UniProtKB/Swiss-Prot:Q02088\] \[protein=tropomyosin\] \[protein_id=CAA93291.2\] \[location=complement(join(4319660..4320127,4320187..4320204))\] \[gbkey=CDS\] | 0.4054527 | 47.17736 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAB59884.1_1634 \[gene=rpp203\] \[locus_tag=SPOM_SPAC1071.08\] \[db_xref=EnsemblGenomes-Gn:SPAC1071.08,EnsemblGenomes-Tr:SPAC1071.08.1,GOA:O14317,InterPro:IPR001859,InterPro:IPR027534,InterPro:IPR038716,PomBase:SPAC1071.08,PomBase:SPAC1071.08.1,UniProtKB/Swiss-Prot:O14317\] \[protein=60S acidic ribosomal protein A2\] \[protein_id=CAB59884.1\] \[location=join(3869145..3869453,3869548..3869571)\] \[gbkey=CDS\] | 0.4042300 | 43.70731 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAK9837800.1_831 \[gene=mzt1\] \[locus_tag=SPOM_SPAC9G1.15C\] \[protein=mitotic spindle organizing protein Mzt1\] \[protein_id=CAK9837800.1\] \[location=complement(1985350..1985544)\] \[gbkey=CDS\] | 0.4027027 | 48.93427 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAA97345.1_1886 \[gene=dhm2\] \[locus_tag=SPOM_SPAC17C9.15C\] \[db_xref=EnsemblGenomes-Gn:SPAC17C9.15c,EnsemblGenomes-Tr:SPAC17C9.15c.1,GOA:Q10486,PomBase:SPAC17C9.15c,PomBase:SPAC17C9.15c.1,UniProtKB/Swiss-Prot:Q10486\] \[protein=Schizosaccharomyces specific protein Dhm2\] \[protein_id=CAA97345.1\] \[location=4475290..4475574\] \[gbkey=CDS\] | 0.4020919 | 48.23675 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAA22608.1_373 \[locus_tag=SPOM_SPAC1687.14C\] \[db_xref=EnsemblGenomes-Gn:SPAC1687.14c,EnsemblGenomes-Tr:SPAC1687.14c.1,GOA:O94455,InterPro:IPR002048,InterPro:IPR011992,PomBase:SPAC1687.14c,PomBase:SPAC1687.14c.1,UniProtKB/Swiss-Prot:O94455\] \[protein=EF hand family protein, centrin-like\] \[protein_id=CAA22608.1\] \[location=complement(928213..928443)\] \[gbkey=CDS\] | 0.3990692 | 48.93603 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAB90142.1_1151 \[gene=rpp101\] \[locus_tag=SPOM_SPAC644.15\] \[db_xref=GOA:P17476,InterPro:IPR027534,InterPro:IPR038716,UniProtKB/Swiss-Prot:P17476\] \[protein=ribosomal protein P1 Rpp101\] \[protein_id=CAB90142.1\] \[location=join(2700453..2700521,2700658..2700918)\] \[gbkey=CDS\] | 0.3976570 | 43.02944 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAA93895.1_2135 \[gene=rrn10\] \[locus_tag=SPOM_SPAC22E12.08\] \[db_xref=EnsemblGenomes-Gn:SPAC22E12.08,EnsemblGenomes-Tr:SPAC22E12.08.1,GOA:Q10360,InterPro:IPR022793,PomBase:SPAC22E12.08,PomBase:SPAC22E12.08.1,UniProtKB/Swiss-Prot:Q10360\] \[protein=RNA polymerase I upstream activation factor complex subunit Rrn10\] \[protein_id=CAA93895.1\] \[location=join(5033478..5033603,5033656..5033823)\] \[gbkey=CDS\] | 0.3929006 | 47.79863 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAK9838068.1_1099 \[gene=rrp36\] \[locus_tag=SPOM_SPAC823.04\] \[protein=rRNA processing protein Rrp36\] \[protein_id=CAK9838068.1\] \[location=join(2587729..2587897,2587952..2588025,2588025..2588066,2588066..2588620)\] \[gbkey=CDS\] | 0.3923268 | 49.85812 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAA93352.1_1625 \[gene=pcc1\] \[locus_tag=SPOM_SPAC4H3.13\] \[db_xref=EnsemblGenomes-Gn:SPAC4H3.13,EnsemblGenomes-Tr:SPAC4H3.13.1,GOA:Q10220,InterPro:IPR015419,PomBase:SPAC4H3.13,PomBase:SPAC4H3.13.1,UniProtKB/Swiss-Prot:Q10220\] \[protein=EKC/KEOPS complex subunit Pcc1\] \[protein_id=CAA93352.1\] \[location=join(3853081..3853112,3853190..3853424)\] \[gbkey=CDS\] | 0.3911036 | 51.55039 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CBA11500.1_1175 \[locus_tag=SPOM_SPAC6F6.19\] \[db_xref=EnsemblGenomes-Gn:SPAC6F6.19,EnsemblGenomes-Tr:SPAC6F6.19.1,GOA:C6Y4A5,InterPro:IPR025239,InterPro:IPR039249,PomBase:SPAC6F6.19,PomBase:SPAC6F6.19.1,UniProtKB/Swiss-Prot:C6Y4A5\] \[protein=RNA-binding protein, G-patch type, human GPATCH11 ortholog\] \[protein_id=CBA11500.1\] \[location=complement(join(2759155..2759292,2759479..2759709))\] \[gbkey=CDS\] | 0.3908086 | 51.03777 | Unannotated / No Pfam Domain Found |

##### Bottom 10 Regulatory/Slowing Translation Genes (Structural Control)

| Protein ID (NCBI accession) | tAI | ENC | GO_Terms |
|:------------------------------------------------------------|--:|--:|:-----|
| lcl\|CU329670.1_cds_CCD31337.1_1757 \[gene=smp1\] \[locus_tag=SPOM_SPAC25B8.20\] \[db_xref=EnsemblGenomes-Gn:SPAC25B8.20,EnsemblGenomes-Tr:SPAC25B8.20.1,GOA:G2TRN5,PomBase:SPAC25B8.20,PomBase:SPAC25B8.20.1,UniProtKB/Swiss-Prot:G2TRN5\] \[protein=Schizosaccharomyces specific microprotein Smp1\] \[protein_id=CCD31337.1\] \[location=join(4166830..4166968,4167077..4167156)\] \[gbkey=CDS\] | 0.3056423 | 50.67601 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CBA11507.1_2039 \[locus_tag=SPOM_SPAC20G4.09\] \[db_xref=EnsemblGenomes-Gn:SPAC20G4.09,EnsemblGenomes-Tr:SPAC20G4.09.1,GOA:C6Y4B1,PomBase:SPAC20G4.09,PomBase:SPAC20G4.09.1,UniProtKB/Swiss-Prot:C6Y4B1\] \[protein=Schizosaccharomyces pombe specific protein\] \[protein_id=CBA11507.1\] \[location=join(4826230..4826297,4826343..4826349)\] \[gbkey=CDS\] | 0.3141471 | 55.72568 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CBA11502.1_1410 \[gene=iec5\] \[locus_tag=SPOM_SPAPB1E7.14\] \[db_xref=EnsemblGenomes-Gn:SPAPB1E7.14,EnsemblGenomes-Tr:SPAPB1E7.14.1,GOA:B8Y7Y5,PomBase:SPAPB1E7.14,PomBase:SPAPB1E7.14.1,UniProtKB/Swiss-Prot:B8Y7Y5\] \[protein=Ino80 complex subunit Iec5\] \[protein_id=CBA11502.1\] \[location=join(3321979..3322018,3322071..3322331,3322379..3322512)\] \[gbkey=CDS\] | 0.3188328 | 48.83618 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAC37429.2_2083 \[gene=min8\] \[locus_tag=SPOM_SPAPB8E5.10\] \[db_xref=EnsemblGenomes-Gn:SPAPB8E5.10,EnsemblGenomes-Tr:SPAPB8E5.10.1,GOA:Q9C0X5,PomBase:SPAPB8E5.10,PomBase:SPAPB8E5.10.1,UniProtKB/Swiss-Prot:Q9C0X5\] \[protein=mitochondrial mini protein, single membrane pass Min8, meiosis specific splicing in fission yeast\] \[protein_id=CAC37429.2\] \[location=join(4926538..4926610,4926659..4926729,4926795..4926872)\] \[gbkey=CDS\] | 0.3208721 | 50.38650 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAA97362.1_2193 \[gene=mug106\] \[locus_tag=SPOM_SPAC26F1.05\] \[db_xref=EnsemblGenomes-Gn:SPAC26F1.05,EnsemblGenomes-Tr:SPAC26F1.05.1,GOA:Q10493,PomBase:SPAC26F1.05,PomBase:SPAC26F1.05.1,UniProtKB/Swiss-Prot:Q10493\] \[protein=Schizosaccharomyces pombe specific protein Mug106\] \[protein_id=CAA97362.1\] \[location=complement(5174979..5175326)\] \[gbkey=CDS\] | 0.3243809 | 53.64945 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAC48261.1_1118 \[gene=ost4\] \[locus_tag=SPOM_SPAC7D4.15C\] \[db_xref=EnsemblGenomes-Gn:SPAC7D4.15c,EnsemblGenomes-Tr:SPAC7D4.15c.1,GOA:Q96VG2,InterPro:IPR018943,InterPro:IPR036330,PomBase:SPAC7D4.15c,PomBase:SPAC7D4.15c.1,UniProtKB/Swiss-Prot:Q96VG2\] \[protein=oligosaccharyltransferase subunit Ost4\] \[protein_id=CAC48261.1\] \[location=complement(2623940..2624038)\] \[gbkey=CDS\] | 0.3250023 | 49.47809 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAB57336.1_708 \[gene=tim14\] \[locus_tag=SPOM_SPAC824.06\] \[db_xref=EnsemblGenomes-Gn:SPAC824.06,EnsemblGenomes-Tr:SPAC824.06.1,GOA:Q9UT37,InterPro:IPR001623,InterPro:IPR036869,PomBase:SPAC824.06,PomBase:SPAC824.06.1,UniProtKB/Swiss-Prot:Q9UT37\] \[protein=TIM23 translocase complex subunit Tim14\] \[protein_id=CAB57336.1\] \[location=join(1695102..1695104,1695152..1695199,1695253..1695624)\] \[gbkey=CDS\] | 0.3264218 | 53.64483 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAB54828.1_2157 \[gene=rpl3002\] \[locus_tag=SPOM_SPAC1250.05\] \[db_xref=EnsemblGenomes-Gn:SPAC1250.05,EnsemblGenomes-Tr:SPAC1250.05.1,GOA:Q9UTP0,InterPro:IPR004038,InterPro:IPR022991,InterPro:IPR029064,InterPro:IPR039109,PomBase:SPAC1250.05,PomBase:SPAC1250.05.1,UniProtKB/Swiss-Prot:Q9UTP0\] \[protein=60S ribosomal protein L30\] \[protein_id=CAB54828.1\] \[location=complement(join(5095211..5095301,5095548..5095810))\] \[gbkey=CDS\] | 0.3270606 | 51.24224 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAB10854.1_622 \[gene=alg14\] \[locus_tag=SPOM_SPAC5D6.06C\] \[db_xref=EnsemblGenomes-Gn:SPAC5D6.06c,EnsemblGenomes-Tr:SPAC5D6.06c.1,GOA:O14199,InterPro:IPR013969,PomBase:SPAC5D6.06c,PomBase:SPAC5D6.06c.1,UniProtKB/Swiss-Prot:O14199\] \[protein=UDP-GlcNAc transferase associated protein Alg14\] \[protein_id=CAB10854.1\] \[location=join(1502909..1503373,1503420..1503587)\] \[gbkey=CDS\] | 0.3277142 | 47.11789 | Unannotated / No Pfam Domain Found |
| lcl\|CU329670.1_cds_CAB11270.1_1307 \[gene=mrp10\] \[locus_tag=SPOM_SPAC24C9.13C\] \[db_xref=EnsemblGenomes-Gn:SPAC24C9.13c,EnsemblGenomes-Tr:SPAC24C9.13c.1,GOA:O13973,InterPro:IPR010625,InterPro:IPR017264,PomBase:SPAC24C9.13c,PomBase:SPAC24C9.13c.1,UniProtKB/Swiss-Prot:O13973\] \[protein=mitochondrial ribosomal protein subunit S37\] \[protein_id=CAB11270.1\] \[location=complement(join(3071685..3071750,3071804..3072004))\] \[gbkey=CDS\] | 0.3289369 | 51.25867 | Unannotated / No Pfam Domain Found |

#### Local Protein Feature Correlation

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-organism-deep-dive-loop-6.png)
\#### Overlap of Structural Domains and Sequence Complexity

| Protein ID (NCBI accession) | Local_Cohort | Overlapping_Domains |
|:-------------------------------------------------------------|:----|:----|
| lcl\|CU329670.1_cds_CAC05742.1_4 \[locus_tag=SPOM_SPAC212.08C\] \[db_xref=EnsemblGenomes-Gn:SPAC212.08c,EnsemblGenomes-Tr:SPAC212.08c.1,GOA:Q9HGP7,PomBase:SPAC212.08c,PomBase:SPAC212.08c.1,UniProtKB/Swiss-Prot:Q9HGP7\] \[protein=S. pombe specific GPI anchored protein family 1\] \[protein_id=CAC05742.1\] \[location=12158..12994\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAC05738.1_9 \[locus_tag=SPOM_SPAC212.04C\] \[db_xref=EnsemblGenomes-Gn:SPAC212.04c,EnsemblGenomes-Tr:SPAC212.04c.1,GOA:Q9HGP8,InterPro:IPR009340,PomBase:SPAC212.04c,PomBase:SPAC212.04c.1,UniProtKB/Swiss-Prot:Q9HGP8\] \[protein=S. pombe specific DUF999 family protein 1\] \[protein_id=CAC05738.1\] \[location=join(21587..22076,22132..22508)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAC05737.1_10 \[locus_tag=SPOM_SPAC212.03\] \[db_xref=EnsemblGenomes-Gn:SPAC212.03,EnsemblGenomes-Tr:SPAC212.03.1,GOA:Q9HGP9,PomBase:SPAC212.03,PomBase:SPAC212.03.1,UniProtKB/Swiss-Prot:Q9HGP9\] \[protein=hypothetical protein\] \[protein_id=CAC05737.1\] \[location=complement(23589..23978)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAC05736.1_11 \[locus_tag=SPOM_SPAC212.02\] \[db_xref=EnsemblGenomes-Gn:SPAC212.02,EnsemblGenomes-Tr:SPAC212.02.1,GOA:Q9HGQ0,PomBase:SPAC212.02,PomBase:SPAC212.02.1,UniProtKB/Swiss-Prot:Q9HGQ0\] \[protein=Schizosaccharomyces pombe specific protein\] \[protein_id=CAC05736.1\] \[location=complement(27353..27763)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAB69633.1_24 \[gene=fex1\] \[locus_tag=SPOM_SPAC977.11\] \[db_xref=GOA:P0CU19,InterPro:IPR003691,UniProtKB/Swiss-Prot:P0CU19\] \[protein=plasma membrane fluoride export channel Fex1\] \[protein_id=CAB69633.1\] \[location=55274..56209\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAB69634.1_25 \[locus_tag=SPOM_SPAC977.12\] \[db_xref=EnsemblGenomes-Gn:SPAC977.12,EnsemblGenomes-Tr:SPAC977.12.1,GOA:Q9UTS7,InterPro:IPR004550,InterPro:IPR006034,InterPro:IPR027473,InterPro:IPR027474,InterPro:IPR036152,InterPro:IPR037152,InterPro:IPR040919,PomBase:SPAC977.12,PomBase:SPAC977.12.1,UniProtKB/Swiss-Prot:Q9UTS7\] \[protein=L-asparaginase\] \[protein_id=CAB69634.1\] \[location=56492..57562\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CCD31308.1_31 \[locus_tag=SPOM_SPAPJ695.02\] \[db_xref=EnsemblGenomes-Gn:SPAPJ695.02,EnsemblGenomes-Tr:SPAPJ695.02.1,GOA:G2TRM2,PomBase:SPAPJ695.02,PomBase:SPAPJ695.02.1,UniProtKB/Swiss-Prot:G2TRM2\] \[protein=Schizosaccharomyces pombe specific protein\] \[protein_id=CCD31308.1\] \[location=complement(75698..75910)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAB03596.1_34 \[gene=shu1\] \[locus_tag=SPOM_SPAC1F8.02C\] \[db_xref=EnsemblGenomes-Gn:SPAC1F8.02c,EnsemblGenomes-Tr:SPAC1F8.02c.1,GOA:Q92340,PomBase:SPAC1F8.02c,PomBase:SPAC1F8.02c.1,UniProtKB/Swiss-Prot:Q92340\] \[protein=cell-surface heme acquisition protein Shu1\] \[protein_id=CAB03596.1\] \[location=complement(85598..86278)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAB03599.1_37 \[gene=isp3\] \[locus_tag=SPOM_SPAC1F8.05\] \[db_xref=EnsemblGenomes-Gn:SPAC1F8.05,EnsemblGenomes-Tr:SPAC1F8.05.1,GOA:P40899,PomBase:SPAC1F8.05,PomBase:SPAC1F8.05.1,UniProtKB/Swiss-Prot:P40899\] \[protein=spore wall structural constituent Isp3\] \[protein_id=CAB03599.1\] \[location=96000..96548\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CCD31309.1_41 \[locus_tag=SPOM_SPAC11D3.19\] \[db_xref=EnsemblGenomes-Gn:SPAC11D3.19,EnsemblGenomes-Tr:SPAC11D3.19.1,GOA:G2TRM3,PomBase:SPAC11D3.19,PomBase:SPAC11D3.19.1,UniProtKB/Swiss-Prot:G2TRM3\] \[protein=Schizosaccharomyces pombe specific protein\] \[protein_id=CCD31309.1\] \[location=106965..107198\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAA92302.1_42 \[locus_tag=SPOM_SPAC11D3.01C\] \[db_xref=EnsemblGenomes-Gn:SPAC11D3.01c,EnsemblGenomes-Tr:SPAC11D3.01c.1,GOA:Q10080,InterPro:IPR018824,PomBase:SPAC11D3.01c,PomBase:SPAC11D3.01c.1,UniProtKB/Swiss-Prot:Q10080\] \[protein=Con-6 family conserved fungal protein\] \[protein_id=CAA92302.1\] \[location=complement(108908..109147)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAA92305.1_45 \[locus_tag=SPOM_SPAC11D3.04C\] \[db_xref=EnsemblGenomes-Gn:SPAC11D3.04c,EnsemblGenomes-Tr:SPAC11D3.04c.1,GOA:Q10083,InterPro:IPR032710,InterPro:IPR037401,PomBase:SPAC11D3.04c,PomBase:SPAC11D3.04c.1,UniProtKB/Swiss-Prot:Q10083\] \[protein=polyketide cyclase SnoaL-like domain protein\] \[protein_id=CAA92305.1\] \[location=complement(112917..113309)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAA92306.1_46 \[gene=mfs2\] \[locus_tag=SPOM_SPAC11D3.05\] \[db_xref=EnsemblGenomes-Gn:SPAC11D3.05,EnsemblGenomes-Tr:SPAC11D3.05.1,GOA:Q10084,InterPro:IPR011701,InterPro:IPR020846,InterPro:IPR036259,PomBase:SPAC11D3.05,PomBase:SPAC11D3.05.1,UniProtKB/Swiss-Prot:Q10084\] \[protein=transmembrane transporter Mfs2\] \[protein_id=CAA92306.1\] \[location=114080..115720\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAA92314.1_54 \[gene=hsp3104\] \[locus_tag=SPOM_SPAC11D3.13\] \[db_xref=EnsemblGenomes-Gn:SPAC11D3.13,EnsemblGenomes-Tr:SPAC11D3.13.1,GOA:Q10092,InterPro:IPR002818,InterPro:IPR029062,PomBase:SPAC11D3.13,PomBase:SPAC11D3.13.1,UniProtKB/Swiss-Prot:Q10092\] \[protein=ThiJ domain protein, glyoxylase III\] \[protein_id=CAA92314.1\] \[location=130329..130997\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAA89952.1_61 \[gene=hsp3102\] \[locus_tag=SPOM_SPAC5H10.02C\] \[db_xref=EnsemblGenomes-Gn:SPAC5H10.02c,EnsemblGenomes-Tr:SPAC5H10.02c.1,GOA:Q09675,InterPro:IPR002818,InterPro:IPR029062,PomBase:SPAC5H10.02c,PomBase:SPAC5H10.02c.1,UniProtKB/Swiss-Prot:Q09675\] \[protein=glyoxylase III Hsp3102\] \[protein_id=CAA89952.1\] \[location=complement(148213..148935)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAA89957.1_66 \[locus_tag=SPOM_SPAC5H10.07\] \[db_xref=EnsemblGenomes-Gn:SPAC5H10.07,EnsemblGenomes-Tr:SPAC5H10.07.1,GOA:Q09678,PomBase:SPAC5H10.07,PomBase:SPAC5H10.07.1,UniProtKB/Swiss-Prot:Q09678\] \[protein=Schizosaccharomyces pombe specific protein\] \[protein_id=CAA89957.1\] \[location=159203..159472\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAA91108.1_87 \[gene=rcn1\] \[locus_tag=SPOM_SPAC13G6.15C\] \[db_xref=EnsemblGenomes-Gn:SPAC13G6.15c,EnsemblGenomes-Tr:SPAC13G6.15c.1,GOA:Q09791,InterPro:IPR006931,PomBase:SPAC13G6.15c,PomBase:SPAC13G6.15c.1,UniProtKB/Swiss-Prot:Q09791\] \[protein=serine/threonine protein phosphatase (calcipressin) regulatory subunit Rcn1\] \[protein_id=CAA91108.1\] \[location=complement(202487..202978)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAA91774.1_92 \[gene=mpc2\] \[locus_tag=SPOM_SPAC24B11.09\] \[db_xref=EnsemblGenomes-Gn:SPAC24B11.09,EnsemblGenomes-Tr:SPAC24B11.09.1,GOA:Q09896,InterPro:IPR005336,PomBase:SPAC24B11.09,PomBase:SPAC24B11.09.1,UniProtKB/Swiss-Prot:Q09896\] \[protein=mitochondrial carrier, pyruvate Mpc2\] \[protein_id=CAA91774.1\] \[location=216483..216839\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAB61651.1_96 \[locus_tag=SPOM_SPAC24B11.14\] \[db_xref=EnsemblGenomes-Gn:SPAC24B11.14,EnsemblGenomes-Tr:SPAC24B11.14.1,PomBase:SPAC24B11.14,PomBase:SPAC24B11.14.1,UniProtKB/Swiss-Prot:Q9UR21\] \[protein=dubious\] \[protein_id=CAB61651.1\] \[location=join(230545..230681,230743..230799,230830..231136)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CU329670.1_cds_CAB55280.1_97 \[gene=hem3\] \[locus_tag=SPOM_SPAC24B11.13\] \[db_xref=EnsemblGenomes-Gn:SPAC24B11.13,EnsemblGenomes-Tr:SPAC24B11.13.1,GOA:Q09899,InterPro:IPR000860,InterPro:IPR022417,InterPro:IPR022418,InterPro:IPR022419,InterPro:IPR036803,PomBase:SPAC24B11.13,PomBase:SPAC24B11.13.1,UniProtKB/Swiss-Prot:Q09899\] \[protein=hydroxymethylbilane synthase Hem3\] \[protein_id=CAB55280.1\] \[location=231655..232665\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |

### Anidulans

#### Intra-Organismal tAI Distribution and Core Outliers

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-organism-deep-dive-loop-7.png)
\##### Top 10 High-Efficiency Expression Target Genes (Elite Pool)

| Protein ID (NCBI accession) | tAI | ENC | GO_Terms |
|:-----------------------------------------------------|---:|---:|:----------|
| lcl\|CH236927.1_cds_EAA57619.1_190 \[protein=hypothetical protein\] \[protein_id=EAA57619.1\] \[location=join(547972..548011,548075..548174,548241..548587,548638..549074)\] \[gbkey=CDS\] | 0.2755074 | 49.33653 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA58161.1_146 \[protein=hypothetical protein\] \[protein_id=EAA58161.1\] \[location=join(415239..415299,415470..415553,415666..415727)\] \[gbkey=CDS\] | 0.2714371 | 47.98608 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA57622.1_193 \[protein=hypothetical protein\] \[protein_id=EAA57622.1\] \[location=complement(join(554847..555121,555450..555482,555540..555741))\] \[gbkey=CDS\] | 0.2603454 | 41.71824 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA58186.1_171 \[protein=hypothetical protein\] \[protein_id=EAA58186.1\] \[location=join(484192..484277,484330..484534)\] \[gbkey=CDS\] | 0.2522965 | 48.01142 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA57840.1_14 \[protein=hypothetical protein\] \[protein_id=EAA57840.1\] \[location=complement(join(36309..36486,36561..36814,36847..36941,37035..37086))\] \[gbkey=CDS\] | 0.2521414 | 47.99477 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA58518.1_214 \[protein=hypothetical protein\] \[protein_id=EAA58518.1\] \[location=complement(627985..630417)\] \[gbkey=CDS\] | 0.2496188 | 38.46278 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA58160.1_145 \[protein=hypothetical protein\] \[protein_id=EAA58160.1\] \[location=join(413617..413625,413704..414099,414162..414248,414306..414335)\] \[gbkey=CDS\] | 0.2483709 | 40.82623 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA57631.1_202 \[protein=hypothetical protein\] \[protein_id=EAA57631.1\] \[location=join(588714..588828,588886..589352,589412..590047)\] \[gbkey=CDS\] | 0.2480683 | 48.90908 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA58549.1_245 \[protein=hypothetical protein\] \[protein_id=EAA58549.1\] \[location=join(735682..735919,735972..737101)\] \[gbkey=CDS\] | 0.2472040 | 43.59347 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA57617.1_188 \[protein=predicted protein\] \[protein_id=EAA57617.1\] \[location=complement(join(542633..542724,542784..542858,542991..543214,543281..543303))\] \[gbkey=CDS\] | 0.2460636 | 53.77020 | Unannotated / No Pfam Domain Found |

##### Bottom 10 Regulatory/Slowing Translation Genes (Structural Control)

| Protein ID (NCBI accession) | tAI | ENC | GO_Terms |
|:------------------------------------------------------|---:|---:|:---------|
| lcl\|CH236927.1_cds_EAA58226.1_341 \[protein=hypothetical protein\] \[protein_id=EAA58226.1\] \[location=join(1035357..1035651,1035716..1035789)\] \[gbkey=CDS\] | 0.2101294 | 55.48290 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA58113.1_98 \[protein=hypothetical protein\] \[protein_id=EAA58113.1\] \[location=complement(join(268364..268660,268708..268875))\] \[gbkey=CDS\] | 0.2103814 | 53.54538 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA57608.1_480 \[protein=hypothetical protein\] \[protein_id=EAA57608.1\] \[location=complement(1459751..1460821)\] \[gbkey=CDS\] | 0.2108477 | 54.13188 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA58284.1_399 \[protein=hypothetical protein\] \[protein_id=EAA58284.1\] \[location=join(1205133..1205832,1205897..1206132)\] \[gbkey=CDS\] | 0.2123241 | 56.59449 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA58304.1_419 \[protein=hypothetical protein\] \[protein_id=EAA58304.1\] \[location=complement(join(1257552..1257554,1257641..1258414,1258515..1258559))\] \[gbkey=CDS\] | 0.2124391 | 54.52944 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA57827.1_1 \[protein=hypothetical protein\] \[protein_id=EAA57827.1\] \[location=complement(join(7633..8567,8623..8941,8989..9067,9109..9140))\] \[gbkey=CDS\] | 0.2124936 | 57.26442 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA57836.1_10 \[protein=hypothetical protein\] \[protein_id=EAA57836.1\] \[location=complement(join(26399..26436,26571..27545,27604..27933,28087..28156,28203..28509,28636..28820))\] \[gbkey=CDS\] | 0.2147563 | 56.21269 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA57706.1_465 \[protein=predicted protein\] \[protein_id=EAA57706.1\] \[location=complement(join(1400109..1400209,1400287..1400414,1400474..1400638,1400744..1400940))\] \[gbkey=CDS\] | 0.2153639 | 53.24225 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA57604.1_476 \[protein=hypothetical protein\] \[protein_id=EAA57604.1\] \[location=complement(join(1450567..1451018,1451120..1451299,1451534..1451585))\] \[gbkey=CDS\] | 0.2154467 | 56.28539 | Unannotated / No Pfam Domain Found |
| lcl\|CH236927.1_cds_EAA57913.1_87 \[protein=hypothetical protein\] \[protein_id=EAA57913.1\] \[location=join(227676..227814,227904..228521,228723..228841)\] \[gbkey=CDS\] | 0.2158715 | 53.89909 | Unannotated / No Pfam Domain Found |

#### Local Protein Feature Correlation

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-organism-deep-dive-loop-8.png)
\#### Overlap of Structural Domains and Sequence Complexity

| Protein ID (NCBI accession) | Local_Cohort | Overlapping_Domains |
|:----------------------------------------------------------|:-----|:------|
| lcl\|CH236927.1_cds_EAA57827.1_1 \[protein=hypothetical protein\] \[protein_id=EAA57827.1\] \[location=complement(join(7633..8567,8623..8941,8989..9067,9109..9140))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57831.1_5 \[protein=hypothetical protein\] \[protein_id=EAA57831.1\] \[location=complement(join(14724..14978,15116..15189,15514..15712))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57836.1_10 \[protein=hypothetical protein\] \[protein_id=EAA57836.1\] \[location=complement(join(26399..26436,26571..27545,27604..27933,28087..28156,28203..28509,28636..28820))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57837.1_11 \[protein=hypothetical protein\] \[protein_id=EAA57837.1\] \[location=complement(join(29375..29520,29574..29618,29674..29857))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57838.1_12 \[protein=hypothetical protein\] \[protein_id=EAA57838.1\] \[location=join(30234..30378,30434..30652,30764..30772,30851..31140,31199..31304,31350..31575,31626..32024,32322..32333,32483..32779,32893..33010,33077..33166,33213..33338)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57839.1_13 \[protein=hypothetical protein\] \[protein_id=EAA57839.1\] \[location=complement(join(33732..33872,33930..34074,34117..34213,34272..34649,34720..34807,34860..34970,35141..35143))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57840.1_14 \[protein=hypothetical protein\] \[protein_id=EAA57840.1\] \[location=complement(join(36309..36486,36561..36814,36847..36941,37035..37086))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57841.1_15 \[protein=hypothetical protein\] \[protein_id=EAA57841.1\] \[location=complement(join(38470..39132,39187..39562,39614..39633))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57845.1_19 \[protein=hypothetical protein\] \[protein_id=EAA57845.1\] \[location=complement(join(50005..50029,50101..50176,50236..50368,50421..50888,50945..51474,51692..51864,51926..52007,52133..52253))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57856.1_30 \[protein=hypothetical protein\] \[protein_id=EAA57856.1\] \[location=complement(join(82096..82104,82163..82197,82290..82641))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57865.1_39 \[protein=FDH_EMENI Probable formate dehydrogenase (NAD-dependent formate dehydrogenase) (FDH)\] \[protein_id=EAA57865.1\] \[location=complement(join(105991..106933,106993..107168,107223..107237))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57867.1_41 \[protein=hypothetical protein\] \[protein_id=EAA57867.1\] \[location=complement(join(112535..113239,113312..113458))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57870.1_44 \[protein=hypothetical protein\] \[protein_id=EAA57870.1\] \[location=complement(join(120954..121438,121495..122401))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57872.1_46 \[protein=hypothetical protein\] \[protein_id=EAA57872.1\] \[location=join(125381..125524,125585..126559)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57874.1_48 \[protein=hypothetical protein\] \[protein_id=EAA57874.1\] \[location=join(136807..137085,137134..137289,137342..137974)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57879.1_53 \[protein=predicted protein\] \[protein_id=EAA57879.1\] \[location=complement(join(147845..147955,148172..148177))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57886.1_60 \[protein=hypothetical protein\] \[protein_id=EAA57886.1\] \[location=complement(join(162673..162806,162927..163140,163196..163234))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57889.1_63 \[protein=hypothetical protein\] \[protein_id=EAA57889.1\] \[location=complement(join(166032..166975,167040..167136))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57897.1_71 \[protein=hypothetical protein\] \[protein_id=EAA57897.1\] \[location=join(188662..188704,188806..188867,189066..189212,189335..189454)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CH236927.1_cds_EAA57898.1_72 \[protein=hypothetical protein\] \[protein_id=EAA57898.1\] \[location=complement(join(189941..190074,190132..190208,190331..190422,190522..190632))\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |

### Calbicans

#### Intra-Organismal tAI Distribution and Core Outliers

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-organism-deep-dive-loop-9.png)
\##### Top 10 High-Efficiency Expression Target Genes (Elite Pool)

| Protein ID (NCBI accession) | tAI | ENC | GO_Terms |
|:--------------------------------------------------------|---:|--:|:--------|
| lcl\|CP017623.1_cds_AOW25982.1_278 \[gene=RPP1A\] \[locus_tag=CAALFM_C103010WA\] \[db_xref=CGD:CAL0000178222\] \[protein=ribosomal protein P1A\] \[protein_id=AOW25982.1\] \[location=630193..630513\] \[gbkey=CDS\] | 0.3938698 | 38.05549 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26730.1_1026 \[gene=RPL29\] \[locus_tag=CAALFM_C111040WA\] \[db_xref=CGD:CAL0000175799\] \[protein=ribosomal 60S subunit protein L29\] \[protein_id=AOW26730.1\] \[location=2433984..2434175\] \[gbkey=CDS\] | 0.3911849 | 46.78706 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26302.1_598 \[gene=RPS14B\] \[locus_tag=CAALFM_C106450CA\] \[db_xref=CGD:CAL0000196584\] \[protein=ribosomal 40S subunit protein S14B\] \[protein_id=AOW26302.1\] \[location=complement(join(1375843..1376231,1376662..1376671))\] \[gbkey=CDS\] | 0.3869592 | 36.98804 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26915.1_1211 \[gene=RPL14\] \[locus_tag=CAALFM_C113050WA\] \[db_xref=CGD:CAL0000179395\] \[protein=ribosomal 60S subunit protein L14B\] \[protein_id=AOW26915.1\] \[location=join(2839748..2839858,2840213..2840497)\] \[gbkey=CDS\] | 0.3815214 | 39.27798 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW25983.1_279 \[gene=RPL13\] \[locus_tag=CAALFM_C103020CA\] \[db_xref=CGD:CAL0000201310\] \[protein=ribosomal 60S subunit protein L13A\] \[protein_id=AOW25983.1\] \[location=complement(631014..631622)\] \[gbkey=CDS\] | 0.3799512 | 39.42993 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26100.1_396 \[gene=HHT21\] \[locus_tag=CAALFM_C104260WA\] \[db_xref=CGD:CAL0000180612\] \[protein=Hht21p\] \[protein_id=AOW26100.1\] \[location=882165..882575\] \[gbkey=CDS\] | 0.3770345 | 40.11821 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26166.1_462 \[locus_tag=CAALFM_C104940CA\] \[db_xref=CGD:CAL0000186239\] \[protein=hypothetical protein\] \[protein_id=AOW26166.1\] \[location=complement(1023215..1023751)\] \[gbkey=CDS\] | 0.3739331 | 38.28586 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW25830.1_126 \[gene=RPS21B\] \[locus_tag=CAALFM_C101370CA\] \[db_xref=CGD:CAL0000190545\] \[protein=ribosomal 40S subunit protein S21B\] \[protein_id=AOW25830.1\] \[location=complement(join(269372..269611,269989..270012))\] \[gbkey=CDS\] | 0.3723464 | 41.96508 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26092.1_388 \[gene=HTA2\] \[locus_tag=CAALFM_C104170CA\] \[db_xref=CGD:CAL0000200441\] \[protein=histone H2A\] \[protein_id=AOW26092.1\] \[location=complement(870671..871066)\] \[gbkey=CDS\] | 0.3716716 | 40.57264 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26538.1_834 \[gene=COX5\] \[locus_tag=CAALFM_C109030CA\] \[db_xref=CGD:CAL0000179426\] \[protein=cytochrome c oxidase subunit Va\] \[protein_id=AOW26538.1\] \[location=complement(1966282..1966779)\] \[gbkey=CDS\] | 0.3706015 | 37.11754 | Unannotated / No Pfam Domain Found |

##### Bottom 10 Regulatory/Slowing Translation Genes (Structural Control)

| Protein ID (NCBI accession) | tAI | ENC | GO_Terms |
|:------------------------------------------------------|---:|---:|:---------|
| lcl\|CP017623.1_cds_AOW26899.1_1195 \[locus_tag=CAALFM_C112840WA\] \[db_xref=CGD:CAL0000187076\] \[protein=hypothetical protein\] \[protein_id=AOW26899.1\] \[location=2797709..2798062\] \[gbkey=CDS\] | 0.2446960 | 46.01996 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26066.1_362 \[gene=SEN15\] \[locus_tag=CAALFM_C103900WA\] \[db_xref=CGD:CAL0000181168\] \[protein=Sen15p\] \[protein_id=AOW26066.1\] \[location=825308..825823\] \[gbkey=CDS\] | 0.2558808 | 41.49212 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26415.1_711 \[gene=HAP2\] \[locus_tag=CAALFM_C107680WA\] \[db_xref=CGD:CAL0000184546\] \[protein=Hap2p\] \[protein_id=AOW26415.1\] \[location=1668121..1669215\] \[gbkey=CDS\] | 0.2609795 | 42.99081 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26146.1_442 \[locus_tag=CAALFM_C104720WA\] \[db_xref=CGD:CAL0000194999\] \[protein=hypothetical protein\] \[protein_id=AOW26146.1\] \[location=982170..982391\] \[gbkey=CDS\] | 0.2622977 | 53.29960 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26234.1_530 \[locus_tag=CAALFM_C105690CA\] \[db_xref=CGD:CAL0000190853\] \[protein=hypothetical protein\] \[protein_id=AOW26234.1\] \[location=complement(1190242..1190988)\] \[gbkey=CDS\] | 0.2628423 | 52.29582 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW25915.1_211 \[locus_tag=CAALFM_C102300WA\] \[db_xref=CGD:CAL0000200924\] \[protein=prefolding complex chaperone subunit\] \[protein_id=AOW25915.1\] \[location=483698..484063\] \[gbkey=CDS\] | 0.2635406 | 41.23686 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26884.1_1180 \[locus_tag=CAALFM_C112690CA\] \[db_xref=CGD:CAL0000191895\] \[protein=hypothetical protein\] \[protein_id=AOW26884.1\] \[location=complement(2766270..2766920)\] \[gbkey=CDS\] | 0.2635439 | 46.69956 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26129.1_425 \[gene=ELC1\] \[locus_tag=CAALFM_C104550WA\] \[db_xref=CGD:CAL0000177304\] \[protein=elongin C\] \[protein_id=AOW26129.1\] \[location=939165..939467\] \[gbkey=CDS\] | 0.2642240 | 47.44146 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW25851.1_147 \[locus_tag=CAALFM_C101610CA\] \[db_xref=CGD:CAL0000179833\] \[protein=hypothetical protein\] \[protein_id=AOW25851.1\] \[location=complement(join(324887..325739,325918..325940))\] \[gbkey=CDS\] | 0.2654277 | 45.01466 | Unannotated / No Pfam Domain Found |
| lcl\|CP017623.1_cds_AOW26438.1_734 \[locus_tag=CAALFM_C107980CA\] \[db_xref=CGD:CAL0000183398\] \[protein=hypothetical protein\] \[protein_id=AOW26438.1\] \[location=complement(1737039..1737770)\] \[gbkey=CDS\] | 0.2654955 | 47.37155 | Unannotated / No Pfam Domain Found |

#### Local Protein Feature Correlation

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-organism-deep-dive-loop-10.png)
\#### Overlap of Structural Domains and Sequence Complexity

| Protein ID (NCBI accession) | Local_Cohort | Overlapping_Domains |
|:---------------------------------------------------------|:------|:-------|
| lcl\|CP017623.1_cds_AOW25705.1_1 \[locus_tag=CAALFM_C100020CA\] \[db_xref=CGD:CAL0000184345\] \[protein=hypothetical protein\] \[protein_id=AOW25705.1\] \[location=complement(4409..4720)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25712.1_8 \[locus_tag=CAALFM_C100100CA\] \[db_xref=CGD:CAL0000197186\] \[protein=cardiolipin synthase\] \[protein_id=AOW25712.1\] \[location=complement(16600..17328)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25715.1_11 \[gene=VPS53\] \[locus_tag=CAALFM_C100130CA\] \[db_xref=CGD:CAL0000194290\] \[protein=Vps53p\] \[protein_id=AOW25715.1\] \[location=complement(19615..21714)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25718.1_14 \[locus_tag=CAALFM_C100160CA\] \[db_xref=CGD:CAL0000200313\] \[protein=hypothetical protein\] \[protein_id=AOW25718.1\] \[location=complement(28286..29488)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25720.1_16 \[gene=RPL16A\] \[locus_tag=CAALFM_C100180WA\] \[db_xref=CGD:CAL0000175593\] \[protein=ribosomal 60S subunit protein L16A\] \[protein_id=AOW25720.1\] \[location=31842..32444\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25721.1_17 \[locus_tag=CAALFM_C100190CA\] \[db_xref=CGD:CAL0000186780\] \[protein=hypothetical protein\] \[protein_id=AOW25721.1\] \[location=complement(32680..33318)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25722.1_18 \[locus_tag=CAALFM_C100200CA\] \[db_xref=CGD:CAL0000183783\] \[protein=hypothetical protein\] \[protein_id=AOW25722.1\] \[location=complement(34099..34446)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25724.1_20 \[gene=PHR2\] \[locus_tag=CAALFM_C100220WA\] \[db_xref=CGD:CAL0000181783\] \[protein=1,3-beta-glucanosyltransferase\] \[protein_id=AOW25724.1\] \[location=37851..39485\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25725.1_21 \[gene=BFA1\] \[locus_tag=CAALFM_C100230CA\] \[db_xref=CGD:CAL0000186201\] \[protein=Bfa1p\] \[protein_id=AOW25725.1\] \[location=complement(40484..41770)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25730.1_26 \[locus_tag=CAALFM_C100320WA\] \[db_xref=CGD:CAL0000193936\] \[protein=retromer subunit\] \[protein_id=AOW25730.1\] \[location=48513..49298\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25736.1_32 \[locus_tag=CAALFM_C100380CA\] \[db_xref=CGD:CAL0000175749\] \[protein=hypothetical protein\] \[protein_id=AOW25736.1\] \[location=complement(54837..57074)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25750.1_46 \[locus_tag=CAALFM_C100520WA\] \[db_xref=CGD:CAL0000184166\] \[protein=phosphoglycerate mutase\] \[protein_id=AOW25750.1\] \[location=81233..81901\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25757.1_53 \[gene=TUF1\] \[locus_tag=CAALFM_C100590WA\] \[db_xref=CGD:CAL0000198818\] \[protein=translation elongation factor Tu\] \[protein_id=AOW25757.1\] \[location=95214..96494\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25765.1_61 \[gene=UGA32\] \[locus_tag=CAALFM_C100670CA\] \[db_xref=CGD:CAL0000195594\] \[protein=Uga32p\] \[protein_id=AOW25765.1\] \[location=complement(114127..115878)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25769.1_65 \[gene=TUB2\] \[locus_tag=CAALFM_C100710CA\] \[db_xref=CGD:CAL0000176403\] \[protein=beta-tubulin\] \[protein_id=AOW25769.1\] \[location=complement(join(123869..125170,125335..125370,125586..125597))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25776.1_72 \[locus_tag=CAALFM_C100790WA\] \[db_xref=CGD:CAL0000193312\] \[protein=hypothetical protein\] \[protein_id=AOW25776.1\] \[location=150607..151707\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25779.1_75 \[locus_tag=CAALFM_C100820WA\] \[db_xref=CGD:CAL0000180272\] \[protein=hypothetical protein\] \[protein_id=AOW25779.1\] \[location=156643..157689\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25788.1_84 \[locus_tag=CAALFM_C100910WA\] \[db_xref=CGD:CAL0000179780\] \[protein=hypothetical protein\] \[protein_id=AOW25788.1\] \[location=join(179823..179870,179946..180041,180129..180558,180645..180829)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25800.1_96 \[locus_tag=CAALFM_C101040WA\] \[db_xref=CGD:CAL0000183037\] \[protein=hypothetical protein\] \[protein_id=AOW25800.1\] \[location=206477..206920\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25802.1_98 \[gene=MBF1\] \[locus_tag=CAALFM_C101060WA\] \[db_xref=CGD:CAL0000194672\] \[protein=Mbf1p\] \[protein_id=AOW25802.1\] \[location=208285..208740\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |

### Intra-Organismal tAI Distribution and Core Outliers

#### Goal

Characterize localized cell investment across individual gene tracks and
identify elite, highly optimized expression targets.

#### Methodology

Render a continuous density curve mapping the internal gene-by-gene tAI
distribution profile of a chosen organism. Isolate the absolute extreme
top 10 and bottom 10 genes, printing them out in a tabular index
accompanied by direct GO descriptions.

#### Critical Observations

Review the distribution modality. Multi-modal distribution curves reveal
clear functional gene stratification, whereas highly skewed
distributions indicate extreme specialization towards housekeeping
genes.

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/single-organism-tai-distribution-1.png)
\#### Top 10 High-Efficiency Expression Target Genes (Elite Pool)

| Protein ID (NCBI accession) | tAI | ENC | GO terms based on Pfam domains obtained by mapping on pfam2go, separated with the pipe symbol ‘\|’ |
|:-----------------------------------------------|--:|--:|:------------------|
| lcl\|CM001196.1_cds_EGP91094.1_572 \[locus_tag=MYCGRDRAFT_78335\] \[db_xref=JGIDB:Mycgr3_78335\] \[protein=hypothetical protein\] \[protein_id=EGP91094.1\] \[location=1906441..1906695\] \[gbkey=CDS\] | 0.8695396 | 50.04355 | NA |
| lcl\|CM001196.1_cds_EGP91548.1_1446 \[locus_tag=MYCGRDRAFT_102767\] \[db_xref=JGIDB:Mycgr3_102767\] \[protein=hypothetical protein\] \[protein_id=EGP91548.1\] \[location=4368971..4369159\] \[gbkey=CDS\] | 0.8681674 | 47.83466 | NA |
| lcl\|CM001196.1_cds_EGP91714.1_1779 \[locus_tag=MYCGRDRAFT_90143\] \[db_xref=JGIDB:Mycgr3_90143\] \[protein=hypothetical protein\] \[protein_id=EGP91714.1\] \[location=5396657..5396980\] \[gbkey=CDS\] | 0.8676202 | 49.40411 | NA |
| lcl\|CM001196.1_cds_EGP91119.1_618 \[locus_tag=MYCGRDRAFT_102256\] \[db_xref=JGIDB:Mycgr3_102256\] \[protein=hypothetical protein\] \[protein_id=EGP91119.1\] \[location=2030498..2030665\] \[gbkey=CDS\] | 0.8673825 | 45.61720 | NA |
| lcl\|CM001196.1_cds_EGP91273.1_913 \[locus_tag=MYCGRDRAFT_83763\] \[db_xref=JGIDB:Mycgr3_83763\] \[protein=hypothetical protein\] \[protein_id=EGP91273.1\] \[location=2856073..2856261\] \[gbkey=CDS\] | 0.8659387 | 52.94876 | NA |
| lcl\|CM001196.1_cds_EGP92619.1_356 \[locus_tag=MYCGRDRAFT_65284\] \[db_xref=JGIDB:Mycgr3_65284\] \[protein=hypothetical protein\] \[protein_id=EGP92619.1\] \[location=complement(join(1324674..1324948,1325030..1325366))\] \[gbkey=CDS\] | 0.8648925 | 53.19686 | NA |
| lcl\|CM001196.1_cds_EGP91636.1_1626 \[locus_tag=MYCGRDRAFT_102855\] \[db_xref=JGIDB:Mycgr3_102855\] \[protein=hypothetical protein\] \[protein_id=EGP91636.1\] \[location=4892793..4893053\] \[gbkey=CDS\] | 0.8643951 | 52.84462 | NA |
| lcl\|CM001196.1_cds_EGP92169.1_1264 \[locus_tag=MYCGRDRAFT_66663\] \[db_xref=InterPro:IPR007109,JGIDB:Mycgr3_66663\] \[protein=hypothetical protein\] \[frame=2\] \[partial=5’\] \[protein_id=EGP92169.1\] \[location=complement(3861310..\>3862198)\] \[gbkey=CDS\] | 0.8640481 | 51.02636 | NA |
| lcl\|CM001196.1_cds_EGP91550.1_1448 \[locus_tag=MYCGRDRAFT_107522\] \[db_xref=JGIDB:Mycgr3_107522\] \[protein=hypothetical protein\] \[protein_id=EGP91550.1\] \[location=join(4374013..4374077,4374181..4374294,4374778..4375100,4375157..4375410,4375460..4375798)\] \[gbkey=CDS\] | 0.8637892 | 49.20338 | NA |
| lcl\|CM001196.1_cds_EGP91429.1_1229 \[locus_tag=MYCGRDRAFT_102643\] \[db_xref=JGIDB:Mycgr3_102643\] \[protein=hypothetical protein\] \[protein_id=EGP91429.1\] \[location=3739650..3739931\] \[gbkey=CDS\] | 0.8636675 | 52.06929 | NA |

#### Bottom 10 Regulatory/Slowing Translation Genes (Structural Control/Folding Delays)

| Protein ID (NCBI accession) | tAI | ENC | GO terms based on Pfam domains obtained by mapping on pfam2go, separated with the pipe symbol ‘\|’ |
|:----------------------------------------------|--:|--:|:------------------|
| lcl\|CM001196.1_cds_EGP91819.1_1984 \[locus_tag=MYCGRDRAFT_79161\] \[db_xref=JGIDB:Mycgr3_79161\] \[protein=hypothetical protein\] \[protein_id=EGP91819.1\] \[location=complement(join(6006196..6006264,6006336..6006407,6006494..6006559))\] \[gbkey=CDS\] | 0.8234634 | 49.13446 | NA |
| lcl\|CM001196.1_cds_EGP92200.1_1195 \[locus_tag=MYCGRDRAFT_89647\] \[db_xref=JGIDB:Mycgr3_89647\] \[protein=hypothetical protein\] \[protein_id=EGP92200.1\] \[location=complement(join(3644229..3644323,3644388..3644474,3644596..3644737,3644793..3644798))\] \[gbkey=CDS\] | 0.8268069 | 51.46783 | NA |
| lcl\|CM001196.1_cds_EGP91658.1_1661 \[locus_tag=MYCGRDRAFT_84191\] \[db_xref=JGIDB:Mycgr3_84191\] \[protein=hypothetical protein\] \[protein_id=EGP91658.1\] \[location=5026169..5026330\] \[gbkey=CDS\] | 0.8294889 | 53.10237 | NA |
| lcl\|CM001196.1_cds_EGP91577.1_1501 \[locus_tag=MYCGRDRAFT_102797\] \[db_xref=InterPro:IPR006461,JGIDB:Mycgr3_102797\] \[protein=hypothetical protein\] \[protein_id=EGP91577.1\] \[location=join(4545681..4545764,4545818..4546252)\] \[gbkey=CDS\] | 0.8301620 | 47.33881 | NA |
| lcl\|CM001196.1_cds_EGP91953.1_1716 \[locus_tag=MYCGRDRAFT_90092\] \[db_xref=JGIDB:Mycgr3_90092\] \[protein=hypothetical protein\] \[protein_id=EGP91953.1\] \[location=complement(5235360..5235686)\] \[gbkey=CDS\] | 0.8314079 | 48.61781 | NA |
| lcl\|CM001196.1_cds_EGP91585.1_1519 \[locus_tag=MYCGRDRAFT_53984\] \[db_xref=InterPro:IPR003213,JGIDB:Mycgr3_53984\] \[protein=hypothetical protein\] \[protein_id=EGP91585.1\] \[location=join(4590613..4590661,4590821..4590936,4590993..4591091)\] \[gbkey=CDS\] | 0.8318188 | 45.80838 | NA |
| lcl\|CM001196.1_cds_EGP91727.1_1804 \[locus_tag=MYCGRDRAFT_102916\] \[db_xref=JGIDB:Mycgr3_102916\] \[protein=hypothetical protein\] \[protein_id=EGP91727.1\] \[location=5459538..5459837\] \[gbkey=CDS\] | 0.8332747 | 51.11220 | NA |
| lcl\|CM001196.1_cds_EGP90871.1_129 \[locus_tag=MYCGRDRAFT_106780\] \[db_xref=JGIDB:Mycgr3_106780\] \[protein=hypothetical protein\] \[protein_id=EGP90871.1\] \[location=577551..578981\] \[gbkey=CDS\] | 0.8345274 | 40.77112 | NA |
| lcl\|CM001196.1_cds_EGP90890.1_170 \[locus_tag=MYCGRDRAFT_102033\] \[db_xref=JGIDB:Mycgr3_102033\] \[protein=hypothetical protein\] \[protein_id=EGP90890.1\] \[location=723354..723716\] \[gbkey=CDS\] | 0.8359453 | 51.99785 | NA |
| lcl\|CM001196.1_cds_EGP90865.1_118 \[locus_tag=MYCGRDRAFT_88698\] \[db_xref=JGIDB:Mycgr3_88698\] \[protein=hypothetical protein\] \[protein_id=EGP90865.1\] \[location=join(536720..536772,536829..537015,537089..537157,537226..537396)\] \[gbkey=CDS\] | 0.8364195 | 46.15146 | NA |

### Local Protein Feature Correlation

#### Goal

Examine localized structural design principles inside a single organism
to check if translational speed optimization matches structural folding
requirements and LCR management.

#### Methodology

Build violin/boxplot charts evaluating protein length, total TM domain
count, total TM domain length, LCR counts, and absolute LCR lengths
against Top 10% vs. Bottom 10% subsets.

#### N/C-Terminal LCR Polarization

Construct bar charts comparing LCR density and specific amino acid
compositions across three distinct functional regions: N-terminus
(0.0-0.25), Internal (0.25-0.75), and C-terminus (0.75-1.0). This
verifies whether low-complexity sequences localized at the N-terminus
correlate with reduced tAI to implement a controlled ribosomal
translation initiation delay.

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/local-protein-feature-correlation-1.png)

### Overlap of Structural Domains and Sequence Complexity

#### Goal

Determine if low complexity sequences (LCRs) directly collide with
functional structural folds within high-efficiency translation
frameworks.

#### Methodology

Calculate the physical intersection length between predicted Pfam domain
bounds and localized LCR coordinates. Flag structural elements
experiencing an intersection overlap of greater than 80% total length
across extreme tAI brackets.

#### Critical Observations

High overlap scores in high-tAI genes point to structural optimization
constraints, whereas high overlap in low-tAI genes could signify
regulatory regions.

| Protein ID (NCBI accession) | Local_Cohort | Overlapping_Domains |
|:---------------------------------------------------------|:------|:-------|
| lcl\|CP017623.1_cds_AOW25705.1_1 \[locus_tag=CAALFM_C100020CA\] \[db_xref=CGD:CAL0000184345\] \[protein=hypothetical protein\] \[protein_id=AOW25705.1\] \[location=complement(4409..4720)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25712.1_8 \[locus_tag=CAALFM_C100100CA\] \[db_xref=CGD:CAL0000197186\] \[protein=cardiolipin synthase\] \[protein_id=AOW25712.1\] \[location=complement(16600..17328)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25715.1_11 \[gene=VPS53\] \[locus_tag=CAALFM_C100130CA\] \[db_xref=CGD:CAL0000194290\] \[protein=Vps53p\] \[protein_id=AOW25715.1\] \[location=complement(19615..21714)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25718.1_14 \[locus_tag=CAALFM_C100160CA\] \[db_xref=CGD:CAL0000200313\] \[protein=hypothetical protein\] \[protein_id=AOW25718.1\] \[location=complement(28286..29488)\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25720.1_16 \[gene=RPL16A\] \[locus_tag=CAALFM_C100180WA\] \[db_xref=CGD:CAL0000175593\] \[protein=ribosomal 60S subunit protein L16A\] \[protein_id=AOW25720.1\] \[location=31842..32444\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25721.1_17 \[locus_tag=CAALFM_C100190CA\] \[db_xref=CGD:CAL0000186780\] \[protein=hypothetical protein\] \[protein_id=AOW25721.1\] \[location=complement(32680..33318)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25722.1_18 \[locus_tag=CAALFM_C100200CA\] \[db_xref=CGD:CAL0000183783\] \[protein=hypothetical protein\] \[protein_id=AOW25722.1\] \[location=complement(34099..34446)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25724.1_20 \[gene=PHR2\] \[locus_tag=CAALFM_C100220WA\] \[db_xref=CGD:CAL0000181783\] \[protein=1,3-beta-glucanosyltransferase\] \[protein_id=AOW25724.1\] \[location=37851..39485\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25725.1_21 \[gene=BFA1\] \[locus_tag=CAALFM_C100230CA\] \[db_xref=CGD:CAL0000186201\] \[protein=Bfa1p\] \[protein_id=AOW25725.1\] \[location=complement(40484..41770)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25730.1_26 \[locus_tag=CAALFM_C100320WA\] \[db_xref=CGD:CAL0000193936\] \[protein=retromer subunit\] \[protein_id=AOW25730.1\] \[location=48513..49298\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25736.1_32 \[locus_tag=CAALFM_C100380CA\] \[db_xref=CGD:CAL0000175749\] \[protein=hypothetical protein\] \[protein_id=AOW25736.1\] \[location=complement(54837..57074)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25750.1_46 \[locus_tag=CAALFM_C100520WA\] \[db_xref=CGD:CAL0000184166\] \[protein=phosphoglycerate mutase\] \[protein_id=AOW25750.1\] \[location=81233..81901\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25757.1_53 \[gene=TUF1\] \[locus_tag=CAALFM_C100590WA\] \[db_xref=CGD:CAL0000198818\] \[protein=translation elongation factor Tu\] \[protein_id=AOW25757.1\] \[location=95214..96494\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25765.1_61 \[gene=UGA32\] \[locus_tag=CAALFM_C100670CA\] \[db_xref=CGD:CAL0000195594\] \[protein=Uga32p\] \[protein_id=AOW25765.1\] \[location=complement(114127..115878)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25769.1_65 \[gene=TUB2\] \[locus_tag=CAALFM_C100710CA\] \[db_xref=CGD:CAL0000176403\] \[protein=beta-tubulin\] \[protein_id=AOW25769.1\] \[location=complement(join(123869..125170,125335..125370,125586..125597))\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25776.1_72 \[locus_tag=CAALFM_C100790WA\] \[db_xref=CGD:CAL0000193312\] \[protein=hypothetical protein\] \[protein_id=AOW25776.1\] \[location=150607..151707\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25779.1_75 \[locus_tag=CAALFM_C100820WA\] \[db_xref=CGD:CAL0000180272\] \[protein=hypothetical protein\] \[protein_id=AOW25779.1\] \[location=156643..157689\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25788.1_84 \[locus_tag=CAALFM_C100910WA\] \[db_xref=CGD:CAL0000179780\] \[protein=hypothetical protein\] \[protein_id=AOW25788.1\] \[location=join(179823..179870,179946..180041,180129..180558,180645..180829)\] \[gbkey=CDS\] | Local Bottom 10% (Slow) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25800.1_96 \[locus_tag=CAALFM_C101040WA\] \[db_xref=CGD:CAL0000183037\] \[protein=hypothetical protein\] \[protein_id=AOW25800.1\] \[location=206477..206920\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |
| lcl\|CP017623.1_cds_AOW25802.1_98 \[gene=MBF1\] \[locus_tag=CAALFM_C101060WA\] \[db_xref=CGD:CAL0000194672\] \[protein=Mbf1p\] \[protein_id=AOW25802.1\] \[location=208285..208740\] \[gbkey=CDS\] | Local Top 10% (Fast) | No Overlap Detected (\>80%) |

Identified Severe Domain-LCR Intersections for Calbicans

------------------------------------------------------------------------

## Molecular Code & Codon Efficiency Space

### Comparative Codon Relative Adaptiveness (*w*<sub>*i*</sub> Distributions)

#### Goal

Decode the underlying organismal translation rules by assessing how
specific codon adaptation weights vary across lifestyles or engineering
groups.

#### Methodology

Compile distribution charts of calculated relative adaptiveness weights
(*w*<sub>*i*</sub>) for all 61 sense codons. Contrast these profiles
across ecological classifications or domains.

#### Critical Observations

Identify shifting optimal codon switches. If distinct biological groups
display contrasting weights for identical codons, it highlights
divergent evolution in their underlying tRNA expression profiles.

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/comparative-codon-weights-1.png)

### Comprehensive Codon Usage Fingerprint (RSCU Heatmap)

#### Goal

Map global codon usage preferences simultaneously across all organisms
to identify evolutionary convergent trends.

#### Methodology

Construct a global two-dimensional matrix containing calculated Relative
Synonymous Codon Usage (RSCU) values across all 61 sense codons for all
samples. Project this matrix using a clustered heatmap with independent
hierarchical tree clustering branches for rows (species) and columns
(codons).

#### Critical Observations

Evaluate whether clustering topology follows strict taxonomic lineages
(phylogeny) or shifts toward ecological lifestyle grouping (niche
adaptation).

![](/home/maxi7524/repositories/tAIpipe/data/tutorial_data/output/aggregated/reports/summary_report_files/figure-markdown_github/rscu-global-heatmap-1.png)

    ## png 
    ##   3
