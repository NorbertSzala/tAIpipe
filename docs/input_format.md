# Input format

`tAIpipe` requires a sample metadata table and input FASTA files for each organism.

## Required sample table

Default path:

```text
config/samples.tsv
```

Required columns:
| column              | description                                                 |
| ------------------- | ----------------------------------------------------------- |
| `sample`            | short unique sample name used in output paths               |
| `species`           | species name                                                |
| `genetic_code`      | NCBI genetic code ID used for CDS translation/codon metrics |
| `domain`            | one of `eukarya`, `bacteria`, `archaea`                     |
| `accession`         | genome accession ID                                         |
| `kingdom`           | broad taxonomy group                                        |
| `phylum`            | phylum-level taxonomy                                       |
| `lifestyle`         | ecological or biological lifestyle category                 |
| `genome`            | genome FASTA filename or glob pattern                       |
| `cds`               | CDS FASTA filename or glob pattern                          |
| `bed`               | BED file with selected regions, if available                |
| `downsampled_fasta` | FASTA file used for fast/testing tRNAscan-SE runs           |
| `first_contig`      | FASTA file containing the first contig                      |
| `include`           | whether the sample should be included                       |

Example:

sample	species	genetic_code	domain	accession	kingdom	phylum	lifestyle	genome	cds	bed	downsampled_fasta	first_contig	include
Spombe	Schizosaccharomyces_pombe	3	eukarya	GCA_000002945.2	Fungi	Ascomycota	nectar_tap_saprotroph	GCA_000002945.2*.fna	GCA_000002945.2*.fna	GCA_000002945.2.bed	GCA_000002945.2_regions.fasta	GCA_000002945.2_first_contig.fna	True

### Notes
genome, cds, downsampled_fasta, and first_contig may use glob patterns.
Each pattern must resolve to exactly one file.
Use include=False to keep a sample in the table but skip it in the workflow.