# Detailed Metrics Reference

## Metrics

### tRNA Adaptation Index (tAI)

#### Biological Purpose: 
Estimates the relative translational elongation speed/efficiency of a gene based on how well its codon configuration matches the cellular tRNA pool availability.

#### Mathematical Concept: 
First, the absolute adaptiveness value ($W_i$) for each codon $i$ is calculated by accounting for all matching tRNAs $j$ (anti-codons) and multiplying their copy numbers by domain-specific wobble-pairing constraints:

$$W_i = \sum_{j=1}^{61} (1 - s_{ij}) G_j$$

Where $G_j$ is the gene copy number of tRNA anticodon $j$, and $s_{ij}$ is the wobble-pairing penalty coefficient (e.g., an exact Watson-Crick match has $s_{ij} = 0$, meaning a pairing efficiency of $1 - 0 = 1$).
  
These absolute values are normalized into relative adaptiveness weights ($w_i$) by dividing by the maximum observed absolute value:
$$w_i = \frac{W_i}{\max(W_i)}$$

The global tAI score for a specific gene $g$ containing $L$ codons is then determined by calculating the geometric mean of its constituent codon weights:
$$\text{tAI}_g = \left( \prod_{k=1}^L w_{i_k} \right)^{\frac{1}{L}}$$

#### Interpretation: 
Scores range continuously between $0$ and $1$. A high tAI score indicates a highly optimized gene that utilizes codons serviced by abundant, highly efficient tRNAs, facilitating rapid translation. A low tAI score points to unoptimized configurations; these can induce local ribosomal stalling or pausing, which is often biologically required to guide precise co-translated protein folding pathways.

***

### Frequency of Optimal Codons (FOP)

#### Biological Purpose: 
Measures the direct proportion of codons inside a gene sequence that are statistically designated as "optimal" or highly efficient for translation.

#### Interpretation: 
Expressed as a fraction between $0$ and $1$. Higher FOP values show that a gene relies extensively on the organism's preferred codon pool to maximize translational throughput.

***

### Codon Adaptation Index (CAI)

#### Biological Purpose: 
Gauges how closely a gene's synonymous codon choices mimic the expression-optimized preferences of a reference set of highly expressed, elite house-keeping genes.

#### Interpretation: 
Values range from $0$ to $1$. Higher values signify strong directional selection towards high-expression setups, a hallmark feature of highly transcribed metabolic or structural genes.

***

### Effective Number of Codons (ENC)

#### Biological Purpose: 
Quantifies the general strength of synonymous codon usage bias (CUB) within a gene sequence, operating completely independently of any external tRNA or expression reference sets.

#### Interpretation: 
Scores map along a scale from $20$ to $61$:

An ENC of **20** signifies extreme, absolute bias, where only one single specific synonymous codon is utilized to code for each amino acid family (all other alternative synonyms are entirely excluded).

An ENC of **61** signifies absolute uniformity, meaning all synonymous codons are utilized equally at random distributions without any preference.

**Rule of Thumb:** Lower ENC values indicate stronger, more pronounced codon usage bias.

***

### Relative Synonymous Codon Usage (RSCU)

#### Biological Purpose: 
Determines if specific individual codons within a synonymous family are overrepresented or avoided relative to a null model of equal, uniform usage among all available alternatives.

#### Mathematical Concept: 
$$\text{RSCU}_{ij} = \frac{X_{ij}}{\frac{1}{n_i} \sum_{j=1}^{n_i} X_{ij}}$$

Where $X_{ij}$ is the raw count of the $j$-th codon for the $i$-th amino acid, and $n_i$ represents the total degeneracy (the number of synonymous options available for that specific amino acid family).

#### Interpretation: 
$\text{RSCU} = 1.0$: Represents the null state of zero bias (the codon is used exactly at its expected random frequency).

$\text{RSCU} > 1.0$: Points to positive selection or overrepresentation of that specific codon.

$\text{RSCU} < 1.0$: Points to underrepresentation or active avoidance of that specific codon.

***

### GC and GC3s Content

#### Biological Purpose: 
Measures overall nucleotide composition alongside focused synonymous third-position composition (`GC3s`).

#### Interpretation: 
While overall GC content tracks global chromosome architectures, **GC3s** focuses exclusively on the third position of synonymous codons. Because third-position mutations are frequently silent (they do not alter the final translated amino acid sequence), GC3s acts as a critical baseline indicator to isolate neutral background mutational pressure and genomic drift from targeted translational selection.