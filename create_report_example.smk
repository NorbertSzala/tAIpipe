# Load configuration from the YAML file
configfile: "config.yaml"

# Helper variables for easier access to paths
DIR_PDB = config["paths"]["data_pdb"]
DIR_SPEC = config["paths"]["data_uniprot_species"]
DIR_KING = config["paths"]["data_uniprot_kingdoms"]
DIR_PLOTS = config["paths"]["plots"]
DIR_TABLES = config["paths"]["tables"]

# List of all specific species and kingdoms from config
SPECIES = config["params"]["species"]
KINGDOMS = config["params"]["kingdoms"]

# --- Main Rule ---
rule all:
    input:
        config["paths"]["report_file"]

# --- Data Acquisition ---

rule download_pdb_metadata:
    """
    Downloads PDB index files and sequences. (as defined in Lab 4 instructions)
    """
    output:
        idx = f"{DIR_PDB}/entries.idx",
        entry_type = f"{DIR_PDB}/pdb_entry_type.txt",
        seqres = f"{DIR_PDB}/pdb_seqres.txt"
    shell:
        "python3 scripts/downloads/download_pdb.py --out {DIR_PDB}"

rule download_uniprot_species:
    """
    Downloads the full proteome for a single species from UniProt.
    """
    output:
        fasta = f"{DIR_SPEC}/{{species}}.fasta"
    params:
        name = "{species}"
    shell:
        "python3 scripts/downloads/download_uniprot.py --species {params.name} --out {DIR_SPEC}"

rule download_uniprot_kingdoms:
    """
    Downloads N randomly selected protein sequences for a given kingdom.
    Uses the random_sample_size from config.
    """
    output:
        fasta = f"{DIR_KING}/{{kingdom}}.fasta"
    params:
        name = "{kingdom}",
        sample_size = config["params"]["random_sample_size"]
    shell:
        "python3 scripts/downloads/download_uniprot.py --kingdom {params.name} --limit {params.sample_size} --out {DIR_KING}"

# --- Lab 4 Analysis (PDB Trends & N-Terminus) ---

rule analyze_pdb_trends:
    input:
        idx = f"{DIR_PDB}/entries.idx",
        seqres = f"{DIR_PDB}/pdb_seqres.txt",
        config = "config.yaml"
    output:
        pie = f"{DIR_PLOTS}/pdb_pie_chart.png",
        interactive_html = f"{DIR_PLOTS}/pdb_trends_interactive.html"
    shell:
        "python3 scripts/analysis/pdb_trends.py --config {input.config}"

rule analyze_pdb_trends_gifs:
    """
    Analyzes PDB trends and generates static and animated plots.
    """
    input:
        idx = f"{DIR_PDB}/entries.idx",
        seqres = f"{DIR_PDB}/pdb_seqres.txt",
        config = "config.yaml"
    output:
        abs_gif = f"{DIR_PLOTS}/pdb_trends_absolute.gif",
        pct_gif = f"{DIR_PLOTS}/pdb_trends_percentage.gif"
    shell:
        "python3 scripts/analysis/pdb_trends_gifs.py --config {input.config}"


# rule analyze_n_terminus:
#     input:
#         fasta_files = expand(f"{DIR_SPEC}/{{s}}.fasta", s=SPECIES)
#     output:
#         nterm_table = f"{DIR_TABLES}/n_terminus_stats.csv",
#     shell:
#         "python3 scripts/analysis/protein_stats.py --mode species "
#         "--input_files {input} --config config.yaml "
#         "--out_plots {DIR_PLOTS} --out_tables {DIR_TABLES}"

# --- Rozdzielona Analiza Lab 5 ---

rule analyze_species_stats:
    input: expand(f"{DIR_SPEC}/{{s}}.fasta", s=SPECIES)
    output:
        f"{DIR_PLOTS}/length_comparison_species.png",
        f"{DIR_PLOTS}/aa_content_E_coli_Human_Yeast.png",
        expand(f"{DIR_PLOTS}/histograms/hist_{{s}}.png", s=SPECIES),
        f"{DIR_TABLES}/species_stats.csv",
        f"{DIR_TABLES}/aa_species_stats.csv"
    shell:
        "python3 scripts/analysis/protein_stats.py --mode species "
        "--input_files {input} --config config.yaml "
        "--out_plots {DIR_PLOTS} --out_tables {DIR_TABLES} "
        "--out_hists {DIR_PLOTS}/histograms"


rule analyze_kingdom_stats:
    input: expand(f"{DIR_KING}/{{k}}.fasta", k=KINGDOMS)
    output:
        f"{DIR_PLOTS}/length_comparison_kingdoms.png",
        f"{DIR_PLOTS}/aa_content_kingdoms.png",
        f"{DIR_TABLES}/kingdom_stats.csv",
        f"{DIR_TABLES}/aa_kingdoms_stats.csv"
    shell:
        "python3 scripts/analysis/protein_stats.py --mode kingdoms "
        "--input_files {input} --config config.yaml "
        "--out_plots {DIR_PLOTS} --out_tables {DIR_TABLES} "

rule analyze_compare_pdb_uniprot:
    input:
        pdb = f"{DIR_PDB}/pdb_seqres.txt",
        kingdoms = expand(f"{DIR_KING}/{{k}}.fasta", k=KINGDOMS)
    output:
        f"{DIR_PLOTS}/length_comparison_pdb_vs_uniprot.png",
        f"{DIR_PLOTS}/aa_content_pdb_vs_uniprot.png",
        f"{DIR_TABLES}/comparison_pdb_uniprot.csv",
        f"{DIR_TABLES}/aa_pdb_vs_uniprot.csv"
    shell:
        "python3 scripts/analysis/protein_stats.py --mode compare_pdb "
        "--input_files {input.pdb} {input.kingdoms} --config config.yaml "
        "--out_plots {DIR_PLOTS} --out_tables {DIR_TABLES} " 
        "--out_hists {DIR_PLOTS}"
# --- Generowanie Raportu ---
rule generate_report:
    input:
        f"{DIR_PLOTS}/pdb_trends_interactive.html",
        f"{DIR_PLOTS}/length_comparison_species.png",
        f"{DIR_PLOTS}/length_comparison_kingdoms.png",
        f"{DIR_PLOTS}/aa_content_kingdoms.png",
        # comparison between pdb and uniprot 
        f"{DIR_PLOTS}/length_comparison_pdb_vs_uniprot.png",
        f"{DIR_PLOTS}/aa_content_pdb_vs_uniprot.png",
        f"{DIR_TABLES}/comparison_pdb_uniprot.csv",
        expand(f"{DIR_PLOTS}/histograms/hist_{{s}}.png", s=SPECIES),
        # gifs
        f"{DIR_PLOTS}/pdb_pie_chart.png",
        f"{DIR_PLOTS}/pdb_trends_absolute.gif",
        f"{DIR_PLOTS}/pdb_trends_percentage.gif",
        # aa (which i forgot ...)
        f"{DIR_TABLES}/aa_species_stats.csv",
        f"{DIR_TABLES}/aa_kingdoms_stats.csv",
        f"{DIR_TABLES}/aa_pdb_vs_uniprot.csv"
    output:
        report = config["paths"]["report_file"]
    params:
        species = SPECIES,
        kingdoms = KINGDOMS,
        sample_size = config["params"]["random_sample_size"],
        plots_dir = DIR_PLOTS,
        tables_dir = DIR_TABLES
    shell:
        "python3 scripts/report/create_md_report.py "
        "--output {output.report} "
        "--species '{params.species}' "
        "--kingdoms '{params.kingdoms}' "
        "--sample_size {params.sample_size} "
        "--plots_dir {params.plots_dir} "
        "--tables_dir {params.tables_dir}"