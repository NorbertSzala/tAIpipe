# Creating a test dataset by downsampling genomes and CDS files

To create a lightweight test dataset, run:

```bash
./workflow/scripts/create_test_dataset.sh \
  -C data/CDS/ \
  -G data/genome/ \
  -A "GCA_000219625.1,GCA_000006255.1,GCA_000002945.2,GCA_000149205.2,GCA_000181695.1,GCA_002104945.1,GCA_016906535.1,GCA_905067625.1"
```

Before running this step, the input directories should contain genome and CDS FASTA files whose names start with the corresponding NCBI accession number, for example:

```txt
GCA_000219625.1*.fna
GCA_000006255.1*.fna
GCA_000002945.2*.fna
```

The script uses the accession numbers provided with `-A` to search for matching genome and CDS files in the directories given by `-G` and `-C`.


| Input    | arguments                                                                     |
| -------- | ----------------------------------------------------------------------------- |
| Argument | Description                                                                   |
| -C       | Path to the directory containing CDS FASTA files                              |
| -G       | Path to the directory containing genome FASTA files                           |
| -A       | Comma-separated list of NCBI accession numbers to include in the test dataset |



## What the script does

For each accession number, the script:

- Searches for one matching genome FASTA file in the genome directory.
- Searches for one matching CDS FASTA file in the CDS directory.
- Extracts the first contig from the genome FASTA. 
- Truncates this contig to the first 1,000,000 bp.
- Tries to extract CDS records corresponding to this first contig by matching the contig ID in CDS FASTA headers.
- If no CDS records match the first contig ID, the script uses a fallback and extracts the first 1,000 CDS records.
- Writes the downsampled genome and CDS files to resources/test_data/.

The output directory structure is:

```txt
resources/test_data/
├── CDS/
│   ├── GCA_000219625.1_first_contig_CDS.fna
│   ├── GCA_000006255.1_first_contig_CDS.fna
│   └── ...
└── genome/
    ├── GCA_000219625.1_first_contig_1Mbp.fna
    ├── GCA_000006255.1_first_contig_1Mbp.fna
    └── ...
```

The downsampled genome files contain only the first contig, limited to 1 Mbp. The CDS files contain either CDS records associated with that first contig, if detectable from FASTA headers, or the first 1,000 CDS records as a fallback. This makes the dataset small enough for quick pipeline testing while preserving realistic FASTA structure.


## Creating metadata for the test dataset

The original `Main_Dataset.tsv `contains metadata for all genomes. For the test dataset, only rows corresponding to selected accession numbers are retained.

First, create the metadata output directory:

```bash
mkdir -p resources/test_data/metadata
```

Then create a filtered metadata file:

```bash
{
  head -n 1 resources/metadata/Main_Dataset.tsv

  for accession in \
    GCA_000219625.1 \
    GCA_000006255.1 \
    GCA_000002945.2 \
    GCA_000149205.2 \
    GCA_000181695.1 \
    GCA_002104945.1 \
    GCA_016906535.1 \
    GCA_905067625.1
  do
    grep "$accession" resources/metadata/Main_Dataset.tsv
  done
} > resources/test_data/metadata/test_metadata.tsv
```

This creates:

`resources/test_data/metadata/test_metadata.tsv`

containing only metadata rows for genomes used in the test dataset.