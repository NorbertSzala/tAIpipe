# Metrics and statistical interpretation

## Gene-level codon metrics

### tRNA Adaptation Index (`tAI`)

`tAI` measures how well a gene's codons match the estimated tRNA pool. The pipeline uses tRNA gene copy numbers as an abundance proxy and `cubar::est_trna_weight()` to apply domain-specific codon–anticodon pairing and wobble penalties.

For codon $i$, an absolute weight is conceptually:

```math
W_i = \sum_j (1-s_{ij})G_j
```

where $G_j$ is the copy number of a compatible tRNA and $s_{ij}$ is its pairing penalty. Relative codon weights are normalized, and gene tAI is their geometric mean across the coding sequence.

**Interpretation:** higher values indicate stronger adaptation to the copy-number-derived tRNA pool. They do not directly measure tRNA expression or translation rate.


### Relative Synonymous Codon Usage (`RSCU`)

For codon $i$ encoding an amino acid with $n$ synonymous codons:

```math
RSCU_i = \frac{X_i}{\frac{1}{n}\sum_{j=1}^{n}X_j}
```

- `RSCU = 1`: expected under equal synonymous use;
- `RSCU > 1`: overrepresented;
- `RSCU < 1`: underrepresented.

The pipeline reports:

- `genome_RSCU`: calculated from all valid CDS;
- `reference_RSCU`: calculated from ribosomal reference CDS.


### Codon Adaptation Index (`CAI`)

CAI is calculated with codon weights estimated from the sample-specific cytosolic ribosomal reference CDS set identified by KofamScan.

For a codon, the relative weight is based on its usage relative to the most used synonymous codon in the reference set. Gene CAI is the geometric mean of these weights. In other words, `CAI` for each sequence is calculated as ratio this sequence's `RSCU` to ribosomal reference sequence `RSCU`.

**Interpretation:** higher values indicate stronger similarity to codon usage in the ribosomal reference genes. CAI values are reference- and genetic-code-dependent.


### Effective Number of Codons (`ENC`)

ENC summarizes synonymous codon-use bias. Typical interpretation:

- values near 20: very strong bias;
- values near 61: weak bias or nearly uniform synonymous usage.

### Expected ENC and `delta_ENC`

The expected ENC under the Wright GC3 model is:

```math
ENC_{expected} = 2 + S + \frac{29}{S^2 + (1-S)^2}
```

where $S=GC3s$.

```math
delta_ENC = ENC_observed - ENC_expected
```

- negative `delta_ENC`: stronger codon bias than expected from GC3s alone;
- positive `delta_ENC`: weaker bias than the model expectation;
- values near zero: broadly compatible with the GC3 expectation.

This is a descriptive deviation, not a formal significance test.


### Frequency of Optimal Codons (`FOP`)

FOP is calculated by `cubar::get_fop()` using the active codon table. It represents the fraction of codons classified as optimal by that implementation.

Do not confuse this value with `reference_optimal_codon` in `codon_profiles.tsv`, which explicitly marks codons with the maximum reference-derived CAI weight within a synonymous family.

### GC and GC3s

- `GC`: fraction of G and C nucleotides across the coding sequence.
- `GC3s`: GC fraction at synonymous third-codon positions as defined by the active genetic code.

GC3s is useful for separating compositional pressure from additional codon-bias effects.

### Amino-acid usage

`*_amino_acid_usage.csv` contains per-gene amino-acid composition derived from codon counts. It is an intermediate matrix and is not merged into the canonical gene table.

## Codon-level tRNA fields

| Field                           | Meaning                                  |
| ------------------------------- | ---------------------------------------- |
| `trna_anticodon`                | Compatible anticodon assigned by `cubar` |
| `trna_id`                       | Amino-acid–anticodon key                 |
| `trna_copy_number` / `ac_level` | Copy-number input used for the codon     |
| `trna_absolute_weight` / `W`    | Unnormalized codon adaptiveness          |
| `tRNA_weight` / `w`             | Relative weight used for tAI             |

## Derived within-genome variables

### `tAI_z`

Standard score calculated separately within each genome:

```math
z = \frac{tAI - \overline{tAI}}{SD(tAI)}
```

It allows gene-level models to compare relative tAI position across genomes with different absolute distributions.

### `tAI_percentile`

Within-genome rank scaled to `[0,1]`. Values near `0` identify low-tAI genes and values near `1` identify high-tAI genes.

## Gene-level mixed models

For every configured binary feature, the model is:

```math
tAI_z ~ feature + available covariates + (1 | sample)
```

Default covariates:

```text
log_protein_length_aa
GC3s
```

The feature p-value is obtained from a likelihood-ratio comparison between the full model and a reduced model without the feature. `estimate` is the fixed-effect difference in within-genome tAI standard deviations for feature-present versus feature-absent genes.

`status=singular_fit` means the model fitted but the random-effect structure was estimated at or near a boundary; interpret such results cautiously.

## Genome-level group tests

For each configured genome metric and grouping variable:

- exactly two groups: Wilcoxon rank-sum test;
- more than two groups: Kruskal–Wallis test.

Effects:

- Wilcoxon: median of the second factor level minus median of the first;
- Kruskal–Wallis: epsilon-squared effect size.

The factor-level order determines the sign of the two-group effect and is reported in `comparison`.

## GO enrichment

Genes are ranked by tAI separately within each genome. High and low tails contain:

```text
ceiling(tail_fraction × number of eligible genes)
```

subject to the optional absolute cap and a non-overlap limit of half the eligible genes.

For every GO term and tail, the script builds a 2×2 table within each genome and combines informative strata using a Cochran–Mantel–Haenszel test.

- `common_odds_ratio > 1`: term enriched in the selected tail;
- `common_odds_ratio < 1`: term depleted in the selected tail;
- confidence interval crossing `1`: effect is not clearly separated from no enrichment.

The background is all eligible non-tail genes, including genes without GO annotation. This prevents the analysis from conditioning only on annotated genes.

## Multiple testing

`q_value` is calculated with the configured `statistics.fdr_method` or `go_enrichment.fdr_method`, defaulting to Benjamini–Hochberg (`BH`). GO p-values are adjusted separately for high- and low-tAI tails.

## Main limitations

- tRNA gene copy number is a proxy, not direct tRNA abundance.
- tAI depends on wobble assumptions and genome annotation quality.
- CAI depends on the quality and size of the KofamScan-derived ribosomal reference.
- ENC and GC3s can be unstable for short or compositionally unusual CDS.
- Group tests do not correct for phylogenetic non-independence.
- GO enrichment quality depends on annotation coverage and consistency across genomes.
