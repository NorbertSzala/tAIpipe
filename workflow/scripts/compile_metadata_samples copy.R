#!/usr/bin/env Rscript

# workflow/scripts/compile_metadata_samples.R
# Skrypt agregujący metadane sekwencji CDS, pobierający termy GO online i łączący wyniki z metrykami tAIpipe.

# Wymagane pakiety instalowane w środowisku conda (np. r-essentials, bioconductor-biostrings, bioconductor-biomart)
suppressPackageStartupMessages({
  library(argparse)
  library(Biostrings)
  library(biomaRt)
  library(stringr)
  library(dplyr)
  library(readr)
  library(purrr)
})

# 1. Parsowanie argumentów wejściowych
parser <- ArgumentParser(description = "Kompilacja metadanych próbek i genów dla tAIpipe")
parser$add_argument("--metadata_dataset", required = TRUE, help = "Ścieżka do pliku dataset_test.tsv")
parser$add_argument("--per-genome-dir", required = TRUE, help = "Katalog główny z wynikami per_genome")
parser$add_argument("--cds-dir", required = TRUE, help = "Katalog z wejściowymi plikami CDS (*.fna)")
parser$add_argument("--gcode-dir", required = TRUE, help = "Katalog z kodami genetycznymi (zostawiony dla kompatybilności)")
parser$add_argument("--output", required = TRUE, help = "Ścieżka wyjściowa dla skonsolidowanego pliku samples.tsv")

args <- parser$parse_args()

# Inicjalizacja połączenia z BioMart (Ensembl Fungi jako domyślna baza dla Twoich danych)
message("Inicjalizacja połączenia online z biomaRt...")
mart <- useMart(biomart = "fungi_mart", host = "https://fungi.ensembl.org")

# Funkcja pomocnicza do pobierania GO online na podstawie listy Protein ID lub UniProt ID
fetch_go_terms_online <- function(protein_ids) {
  tryCatch({
    # Usunięcie pustych wartości
    valid_ids <- protein_ids[!is.na(protein_ids) & protein_ids != ""]
    if (length(valid_ids) == 0) return(list())
    
    # Odpytanie Ensembl BioMart
    # Uwaga: atrybuty mogą się różnić w zależności od wybranego datasetu, 
    # dlatego jako uniwersalny fallback stosujemy wyszukiwanie po uniproksb
    ensembl_mart <- useDataset("uniprot", mart = mart) # Przykładowe użycie dedykowanego datasetu
    
    go_data <- getBM(
      attributes = c('protein_id', 'go_id'),
      filters = 'protein_id',
      values = valid_ids,
      mart = ensembl_mart
    )
    
    # Agregacja po protein_id (sklejanie średnikami)
    go_mapped <- go_data %>%
      filter(go_id != "") %>%
      group_by(protein_id) %>%
      summarise(go_list = paste(unique(go_id), collapse = ";"), .groups = "drop")
    
    return(setNames(go_mapped$go_list, go_mapped$protein_id))
  }, error = function(e) {
    message("Ostrzeżenie: Nie udało się pobrać GO przez biomaRt (błąd sieci/brak ID). Użyty zostanie fallback z nagłówków.")
    return(list())
  })
}

# 2. Wczytanie aktywnego arkusza dataset
dataset_df <- read_tsv(args$metadata_dataset, show_col_types = FALSE) %>%
  filter(include == TRUE | include == "True")

if (nrow(dataset_df) == 0) {
  stop("Brak próbek z flagą include=True w pliku dataset.")
}

all_samples_metadata <- list()

# 3. Pętla przetwarzania dla każdej próbki
for (i in 1:nrow(dataset_df)) {
  current_sample <- dataset_df$sample[i]
  cds_pattern    <- dataset_df$cds[i]
  
  message(paste("Przetwarzanie próbki:", current_sample))
  
  # Znajdź plik CDS pasujący do wzorca z dataset.tsv
  cds_files <- Sys.glob(file.path(args$cds_dir, cds_pattern))
  if (length(cds_files) == 0) {
    warning(paste("Brak pliku CDS dla wzorca:", cds_pattern))
    next
  }
  cds_file <- cds_files[1]
  message(paste("Znaleziono plik CDS:", cds_file)) # DODAJ TO
  
  # Odczyt nagłówków FASTA za pomocą Biostrings
  fasta_headers <- names(readDNAStringSet(cds_file))
  
  # Parsowanie skomplikowanych metadanych z nagłówka NCBI za pomocą wyrażeń regularnych
  parsed_cds <- map_df(fasta_headers, function(header) {
    
    # Ekstrakcja kluczowych tagów
    gene_id      <- str_match(header, "\\[locus_tag=([^\\]]+)\\]")[,2]
    protein_id   <- str_match(header, "\\[protein_id=([^\\]]+)\\]")[,2]
    protein_name <- str_match(header, "\\[protein=([^\\]]+)\\]")[,2]
    location     <- str_match(header, "\\[location=([^\\]]+)\\]")[,2]
    # uniprot_id   <- str_match(header, "UniProtKB/Swiss-Prot:([^,\\]]+)")[,2]
    # uniprot_id   <- str_match(header, "UniProtKB(?:/Swiss-Prot|/TrEMBL)?:([A-Za-z0-9]+)")[,2]
    db_xref <- str_match(header, "\\[db_xref=([^\\]]+)\\]")[,2]
    uniprot_id <- str_match(header, "UniProtKB(?:/Swiss-Prot|/TrEMBL)?:([A-Za-z0-9]+)")[,2]
    if(is.na(uniprot_id)) {
      uniprot_id <- str_match(header, "UniProtKB(?:/Swiss-Prot|/TrEMBL)?:([A-Za-z0-9]+)")[,2]
    }


    if(is.na(uniprot_id)) uniprot_id <- str_match(header, "GOA:([^,\\]]+)")[,2] # alternatywne szukanie bazy
    
    # Wyciąganie GO bezpośrednio z nagłówka jako wbudowany fallback (z db_xref)
    # go_fallback <- paste(str_match_all(header, "GO:([0-9]+)")[[1]][,1], collapse = ";")
    go_matches  <- str_match_all(db_xref, "GO:([0-9]+)")[[1]]
    go_terms    <- if(nrow(go_matches) > 0) paste(unique(go_matches[,2]), collapse = ";") else NA_character_

    # Obliczanie długości na podstawie pozycji w sekwencji (np. complement(<1..5662) lub 4409..4720)
    # Wyciągamy wszystkie cyfry z pola location
    coords <- as.numeric(str_extract_all(location, "[0-9]+")[[1]])
    length_nt <- NA
    length_aa <- NA
    
    if (length(coords) >= 2) {
      length_nt <- max(coords) - min(coords) + 1
      length_aa <- floor((length_nt - 3) / 3) # Odjęcie kodonu stop i podział przez 3
    }
    
    # Główny identyfikator w plikach _summary.tsv (z reguły to cały ciąg przed pierwszą spacją)
    seq_id_metrics <- str_split(header, " ")[[1]][1]
    
    tibble(
      seq_id_metrics = seq_id_metrics,
      gene_id      = ifelse(is.na(gene_id), protein_id, gene_id),
      protein_id   = protein_id,
      uniprot_id   = uniprot_id,
      protein_name = protein_name,
      length_nt    = length_nt,
      length_aa    = length_aa,
      go_fallback  = go_terms,
    )
  })

  # TODO-debugging
  message(paste("Liczba sparsowanych genów w", current_sample, ":", nrow(parsed_cds))) # DODAJ TO
  message("Przykładowe ID z parsed_cds:") # DODAJ TO
  print(head(parsed_cds$seq_id_metrics)) # DODAJ TO
  
  # Próba dociągnięcia termów GO online dla całej próbki
  go_online_map <- fetch_go_terms_online(parsed_cds$protein_id)
  
  # Mapowanie końcowe GO (Priorytet: Online -> Fallback z nagłówka)
  parsed_cds <- parsed_cds %>%
    mutate(
      go_online = as.character(go_online_map[protein_id]),
      go_terms  = case_when(
        !is.na(go_online) & go_online != "" ~ go_online,
        go_fallback != "" ~ go_fallback,
        TRUE ~ NA_character_
      )
    ) %>%
    select(-go_online, -go_fallback)
  
  # Wczytanie powiązanego pliku _summary.tsv wygenerowanego w poprzednich krokach tAIpipe
  summary_file <- file.path(args$per_genome_dir, current_sample, "codon_metrics", paste0(current_sample, "_summary.tsv"))
  
  if (file.exists(summary_file)) {
    metrics_df <- read_tsv(summary_file, show_col_types = FALSE)
    
    #TODO - debugging
    message(paste("Wierszy w parsed_cds:", nrow(parsed_cds)))
    message(paste("Wierszy w metrics_df:", nrow(metrics_df)))
    message("Przykładowe ID z metrics_df (summary):")
    print(head(metrics_df$seq_id))

    metrics_df <- metrics_df %>%
      mutate(seq_id = str_trim(str_extract(seq_id, "^[^\\[]+")))

    # Połączenie sparsowanych metadanych z wyliczonymi metrykami (ENC, CAI, tAI itp.)
    sample_consolidated <- parsed_cds %>%
      inner_join(metrics_df, by = c("seq_id_metrics" = "seq_id")) %>%
      mutate(sample = current_sample) %>%
      select(sample, gene_id, protein_id, uniprot_id, protein_name, length_nt, length_aa, go_terms, 
             ENC, CAI, FOP, tAI, GC, GC3s)

    #TODO - debugging:
    message(paste("Liczba wierszy po join:", nrow(sample_consolidated)))
    
    all_samples_metadata[[current_sample]] <- sample_consolidated
  } else {
    warning(paste("Brak pliku podsumowania miar:", summary_file))
  }
}

# Po pętli:
full_df <- bind_rows(all_samples_metadata)

# Znajdź wszystkie brakujące UniProt
missing_ids <- full_df %>% 
  filter(is.na(uniprot_id)) %>% 
  pull(protein_id) %>% unique()

# Jeśli mamy braki, pobierz je masowo z BioMart
if (length(missing_ids) > 0) {
  message("Pobieranie brakujących danych z bazy...")
  # Użyj swojej funkcji fetch_go_terms_online, ale zmodyfikuj ją, 
  # by zwracała też mapowanie protein_id -> uniprot_id
}

# 4. Zapisanie skonsolidowanej tabeli końcowej dla wszystkich uwzględnionych próbek
final_metadata_df <- bind_rows(all_samples_metadata)

ensure_dir <- dirname(args$output)
if (!dir.exists(ensure_dir)) dir.create(ensure_dir, recursive = TRUE)

readr::write_tsv(final_metadata_df, args$output)
message(paste("Sukces! Metadane zapisano pomyślnie do:", args$output))