#```{r libraries}
library(dplyr)
library(tidyr)
library(rlang)
library(ggplot2)
#```

#```{r TO DO -- DONE}
#rozwiazalem problem z duplikowaniem rzedow. Niektore bialka moga miec kilka takich samych domen. Dotychczasowo rozdzielam takie domeny na tyle samo rzedow. Moj kod do tej pory zliczal liczbe takich domen, a nie liczbe bialek. Teraz usune duplikaty z finalnych tabeli i bedzie w porzadku

############ DONE #############
# DONE    jedna tabela skupiajaca 5 lub 10% top/bottom bialek pod kątem wartości tAI. Z niej opisuje ile z nich mialo ile LCR, jakiej długości, ile TM itp. Comment and PFAM dodaj tak jak do pozostałych
#Wynik (ten sam lub inny): spośród 1000 bialek o najwyzszym tAI 38 z nich ma peptydy sygnalowe

#nastepne tabelki maja byc w formie:
#   kazda tabelka dla jednej cechy np. top and bottom 5%/10% dla bialek z obecnością LCR, signal itd.
#   W tych tabelkach zobacz, jakie PFAMy są obecne w bialkach o skrajnych wartościach tAI ktore mają jakąś cechę
#```

#```{r TO DO -- DONE}
# ########### PYTANIE
#  
# Phy	phyGlomero ..						
# 	- Co łączy domeny tAI dependent i te indepenedent (wszędzie)						
# 	Liczebność tych tAIów powyżej 0.6, 0.7						
# 	Zobaczyć jakie domeny / organizmy są w 1% hi low tAI -> porównamy do  tego zbioru						
# 	- Porównać liczebność % z całym procentem ile braliśmy pod uwagę / analizowaliśmy sekwencji z danego Phylum						
# 	- najpierw bezwzględne zliczenia a potem enrichment. Enrichment policzyłem: (liczba wystąpień domen PFAM = PF00004 w zbiorze extremum/ liczba wszystkich domen w zbiorze extremum hi lub lo)/ (liczba wystąpień domen PFAM = PF00004 w ogólnym zbiorze) / 						
# 	# 	
# 1. zaktualizowac tabele o pfamy
# 2. porownac pfamy -> wnioski opisowe
# 3. powtorzyc ta analize dla 1% i zrobic taka sama tabele +policzyc ile jest bialek w 1%
 
#```


#```{r modify data}
#I need to modify data in another way becouse manipulating data by splitting and merging them is not efficiency and make many mistakes. 

#read whole data
# local path
main_path <- '/home/norbert/IBB/skrypty/Main_Dataset_splitted/data/'
output_path<- '/home/norbert/IBB/skrypty/Main_Dataset_splitted/data/hilo_pfam_ten_percent/'
# output_path<- '/home/norbert/IBB/skrypty/Main_Dataset_splitted/data/hilo_pfam_one_percent/'



#server mycelia path
# main_path <- '/home/norbert_s/tAI/all_chunks/Main_Dataset_splitted/data/'
# output_path<-'/home/norbert_s/tAI/all_chunks/Main_Dataset_splitted/plots_extensive/'#


#Read data with PFAM accession ID and description
PFAM_acc_de <- read.csv(paste0(main_path, 'PFAM_acc_de.csv'), header = FALSE, col.names = c('accession', 'description'))


#już tutaj były niektóre dane zduplikowane. Np. prot_id ==KTW27265.1 zawierał PFAMYY: "PF02349.19 PF02349.19 PF02349.19 PF02349.19 PF12373.12". Niekróre białka mogą mieć zduplikowane domeny i jest to ok.
all_values_table <- read.csv(paste0(main_path, 'shorted_main_dataset.tsv'), sep = '\t', header = FALSE, col.names = c('acc_id', 'prot_id','prot_len', 'domain', 'TM_count', 'TM_length', 'signal', 'LCR_number', 'LCR_length', 'PFAM_LCR', 'PFAM', 'GOterms', 'tAI'))


Phylum <- read.csv(paste0(main_path, 'org_name_phylum.tsv'), sep = '\t', header = TRUE)

all_values_table <- left_join(all_values_table, Phylum, by = c('acc_id' = 'assemblyid'))

#make two sets data. First package can have no described PFAM values. Second, what is more important, is with only that proteins that have assigned PFAMs

all_values_table$PFAM <- trimws(all_values_table$PFAM)
 
#czyli tutaj trzeba zmienic by dodawać opisy PFAMów 
#Sprobuje teraz sumować oryginalne bialka zamiast zduplikowanych
#split PFAM column to make each PFAM id in unique row. The same for goterms
all_values_table <- all_values_table %>%
  separate_rows(PFAM, sep = " ")
all_values_table <- unique(all_values_table)
# all_values_table <- distinct(all_values_table)
  # separate_rows(GOterms, sep = '\\|')
#remove version symbol from PFAM ID
all_values_table$PFAM <- sub("\\..*$", "", all_values_table$PFAM)
all_values_table <- left_join(all_values_table, PFAM_acc_de, by = c('PFAM' = 'accession'))


all_values_only_PFAMs <- all_values_table %>%
  filter(!is.na(PFAM) & PFAM != "")




#table with selected the most important PFAM data
PFAM_shorted <- all_values_only_PFAMs[, c('acc_id', 'prot_id', 'PFAM', 'description', 'GOterms', 'tAI', 'Phylum')]
  
#table with data (must contain PFAM and GOterms simultanously)
GOterms_shorted <- PFAM_shorted %>% filter(GOterms != '')
GOterms_shorted <- GOterms_shorted[, c('acc_id', 'prot_id', 'PFAM', 'description', 'GOterms', 'tAI', 'Phylum')]


#table with proteins that have domain
domain_bolean <- all_values_only_PFAMs[, c('acc_id', 'prot_id', 'PFAM', 'description', 'GOterms', 'tAI', 'domain', 'Phylum')]
domain_bolean <- domain_bolean %>%
    filter(domain == 1) %>% 
    mutate(presence = factor(case_when(
        domain == 1 ~ 'Domain presence'
    ), levels = c('Domain presence')))



#table with proteins that have signal peptides
signal_bolean <- all_values_only_PFAMs[, c('acc_id', 'prot_id', 'PFAM', 'description', 'GOterms', 'tAI', 'signal', 'Phylum')]
signal_bolean <- signal_bolean %>%
  filter(signal == 1) %>% 
  mutate(presence = factor(case_when(
    signal == 1 ~ 'Signal presence',
    ), levels=c('Signal presence')))


#table with proteins and their length
protein_length <- all_values_only_PFAMs[, c('acc_id', 'prot_id', 'PFAM', 'description', 'GOterms', 'tAI', 'prot_len', 'Phylum')]
protein_length <- protein_length %>% 
    mutate(prot_len  = as.numeric(prot_len)) %>% 
    arrange(prot_len) %>% 
    mutate(bin = ntile(prot_len, 10))
colnames(protein_length)[colnames(protein_length) == 'prot_len'] <- 'length'



#table with proteins that have LCR + counted
LCR_number <- all_values_only_PFAMs[, c('acc_id', 'prot_id', 'PFAM', 'description', 'GOterms', 'tAI', 'LCR_number', 'Phylum')]
colnames(LCR_number)[colnames(LCR_number)=='LCR_number'] <- 'LCR_count'
LCR_number <- LCR_number %>%
  filter(LCR_count >= 1) %>% 
  mutate(presence = factor(case_when(
    LCR_count == 1 ~ 'LCR presence',
    TRUE ~ 'LCR absence'), levels = c('LCR presence', 'LCR absence'))) %>% 
  mutate(bin = ntile(LCR_count, 10))



#table with proteins that have domain and LCR
PFAM_LCR <- all_values_only_PFAMs[, c('acc_id', 'prot_id', 'PFAM', 'description', 'GOterms', 'tAI', 'PFAM_LCR', 'Phylum')]
colnames(PFAM_LCR)[colnames(PFAM_LCR)=='PFAM_LCR'] <- 'PFAM_LCR_col'
PFAM_LCR <- PFAM_LCR %>% 
  filter(PFAM_LCR_col != '') %>% 
  mutate(presence = factor(case_when(
    PFAM_LCR_col != '' ~ 'PFAM LCR presence',
    TRUE ~ 'PFAM LCR absence'),
    levels=c('PFAM LCR presence', 'PFAM LCR absence'))) %>% 
  separate_rows(PFAM_LCR_col, sep = ';')




#table with proteins that have Transmembrane elements + length counted
TM_length <- all_values_only_PFAMs[, c('acc_id', 'prot_id', 'PFAM', 'description', 'GOterms', 'tAI', 'TM_length', 'Phylum')]
colnames(TM_length)[colnames(TM_length)=='TM_length'] <- 'TM_length_col'
TM_length <- TM_length %>% 
  filter(TM_length_col > 0) %>% 
  mutate(TM_length_col = as.numeric(TM_length_col)) %>% 
  arrange(TM_length_col) %>% 
  mutate(bin = ntile(TM_length_col, 10))




#table with proteins that have Transmembrane elements
TM_count <- all_values_only_PFAMs[, c('acc_id', 'prot_id', 'PFAM', 'description', 'GOterms', 'tAI', 'TM_count', 'Phylum')]
colnames(TM_count)[colnames(TM_count)=='TM_count'] <- 'TM_count_col'
TM_count <- TM_count %>% 
  filter(TM_count_col > 0) %>%
  mutate(presence = factor(case_when(
    TM_count_col != 0 ~ 'Transmembrane elements presence'),
    levels = c('Transmembrane elements presence'))) %>% 
  mutate(bin = ntile(TM_count_col, 10))



#table with proteins that have domain
LCR_length <- all_values_only_PFAMs[, c('acc_id', 'prot_id', 'PFAM', 'description', 'GOterms', 'tAI', 'LCR_length', 'Phylum')]
LCR_length <- LCR_length %>% 
  filter(LCR_length >= 1) %>% 
    mutate(LCR_length = as.numeric(LCR_length)) %>% 
    arrange(LCR_length) %>% 
    mutate(bin = ntile(LCR_length, 10))



mean_tAI <- mean(domain_bolean$tAI, na.rm=TRUE)
mean_tAI <- mean(LCR_length$tAI, na.rm=TRUE)


#```



#```{r functions}
create_regression_equation <- function(data, formula) {
  model <- lm(formula, data=data)
  summary_model <- summary(model)
  r_sq <- round(summary_model$r.squared, 5)
  coef <- coef(model)
  intercept <- round(coef[1], 2)
  slope <- round(coef[2], 2)
  p_value <- round(summary_model$coefficients[2, 4], 5)
  if (p_value >= 0.05) {
      equation_label <- paste("y = ", intercept, " + ", slope, " * x,", "        R² = ", r_sq, ",        p value ≥ 0.05.", sep = "")
  } else {
      equation_label <- paste("y = ", intercept, " + ", slope, " * x,", "        R² = ", r_sq, ",        p value < 0.05. ", sep = "")

  }
  return(equation_label)
}

#function to create label with statistical test result. Test can be chosen by setting argument parametric_assumption to TRUE (anova or kruskal wallis) or FALSE to T-student or Mann-Whitney.
create_statistical_label <- function(data, var1, var2, parametric_assumption = TRUE) {
  
  if (parametric_assumption) {
    ad_test <- ad.test(var1)#Anderson Darling test to check normality distribution
    
    if (ad_test$p.value < 0.05) {
      # ANOVA test for parametric values with normal distribution
      test = 'ANOVA '
      anova_results <- summary(aov(var1 ~ var2, data = data))
      p_value <- anova_results[[1]]$'Pr(>F)'[1]
    } else {
      # Kruskal-Wallis test for parametric values and lack of normal distribution
      test = 'Kruskal-Wallis test '
      kw_test <- kruskal.test(var1 ~ var2, data = data)
      p_value <- kw_test$p.value
    }
  } else {
    ad_test <- ad.test(var1)
    
    if (ad_test$p.value < 0.05) {
      # Mann-Whitney test #non parametrical test with normal distribution
      test = 'Mann-Whitney test '
      mann_whitney_test <- wilcox.test(var1 ~ var2, data = data)
      p_value <- mann_whitney_test$p.value
    } else {
      # T-test for not normal test not parametrical
      test = 'T-student test '
      t_student_result <- t.test(var1, mu = 0.5)
      p_value <- t_student_result$p.value
    }
  }
  
  # Create the equation label
  if (p_value >= 0.05) {
    equation_label <- paste(test, "p-value ≥ 0.05")
  } else {
    equation_label <- paste(test, "p-value < 0.05")
  }
    

  return(equation_label)
}



get_extreme_tAI_values <- function(table, column, extremum, variable) {
  #' Get Extreme tAI Values
  #' This function returns a dataframe with proteins that have extreme tAI values from top bins 
  #' (the lowest and highest) depending on the bins. The 'extremum' parameter 
  #' is an integer representing the top 'extremum' percentage of values.
  #'
  #' @param table A data frame containing the data.
  #' @param column The name of the column containing tAI values.
  #' @param extremum The top 'extremum' percentage of values.
  #' @param variable Describes whether to consider percentage or bins.
  #' @return A data frame with proteins having extreme tAI values.
  #' @export
  
  # Ensure that the column is numeric
  table[[column]] <- as.numeric(table[[column]])
  
  # Create breaks for binning
  breaks <- seq(0, max(table[[column]], na.rm = TRUE), length.out = extremum + 1)
  
  # If we are using bins
  if (variable == 'bins') {
    table$label <- cut(
      table[[column]],
      breaks = breaks,
      labels = 1:extremum,
      right = FALSE
    )
    
    # Filter the extreme values from the first and last bins
    extremum_table <- table %>% 
      filter(label == 1 | label == extremum)
    
    table$label <- as.character(table$label)  # Convert bin to character so we can change values
    table$label[table$label == '1'] <- 'lowest'  # Replace 1 with 'lowest'
    table$label[table$label == as.character(extremum)] <- 'highest'  # Replace extreme bin with 'highest'
    
    # Filter the extreme values (lowest and highest)
    extremum_table <- table %>% 
      filter(label == 'lowest' | label == 'highest')

  }
    
  if (variable == 'percent') {
    # Sort the values in the column
    sorted_column <- sort(table[[column]], na.last = NA)
    
    # Calculate the lower and upper index for the extremum percentage
    lower_index <- ceiling(length(sorted_column) * (extremum / 100)) # smallest 'extremum'% values
    upper_index <- floor(length(sorted_column) * (1 - (extremum / 100))) # largest 'extremum'% values
  
    # Filter the extreme values
    extremum_table <- table %>%
      filter(table[[column]] <= sorted_column[lower_index] | table[[column]] >= sorted_column[upper_index])
    
    # Add 'lowest' for the smallest 'extremum'% values, 'highest' for the largest
    extremum_table$label <- ifelse(
      extremum_table[[column]] <= sorted_column[lower_index], 
      'lowest', 
      'highest'
    )
  }
  
  return(extremum_table)
}


# function that choose X the most common proteins. Returns these proteins with counted occurences of PFAM(FALSE)/GOterms(TRUE)
count_extreme_goterms <- function(table, x, GOterms = TRUE, extreme_column = NULL) {

  # Determine the sorting column
  sorting_value <- if (GOterms) "GOterms" else "PFAM"
  sorting_sym <- sym(sorting_value)  # convert string to symbol for tidy evaluation

  # Step 1: Summarize counts per label and sorting column (label always contains there values 'highest' or 'lowest' )
  goterm_summary <- table %>%
    group_by(label, !!sorting_sym) %>%
    summarise(n_count = n(), .groups = "drop")
  ##### wydaje mi sie ze tutaj jest blad

  # Step 2: Select the top x most frequent entries for each label
  top_goterms <- goterm_summary %>%
    group_by(label) %>%
    slice_max(n_count, n = x, with_ties = FALSE) %>%
    ungroup()

  # Step 3: Join with the original table to get full info
  result <- table %>%
    semi_join(top_goterms, by = c("label", sorting_value)) %>% 
    left_join(top_goterms, by = c("label", sorting_value))

  # Step 4: Dynamically filter only existing columns
  base_columns <- c("PFAM",'prot_id','description', "GOterms", "tAI", 'mean_tAI', "label", "n_count", "bin", 'Phylum')
  existing_columns <- base_columns[base_columns %in% colnames(result)]  # Keep only existing columns

  if (!is.null(extreme_column) && extreme_column %in% colnames(result)) {
    existing_columns <- c(existing_columns, extreme_column)
  }

  result <- result %>%
    select(all_of(existing_columns))

  return(result)

}



#how to interpret results:
#from certein PFAM with its description there is X proteins (X = n_count for all mentioned bins). It was intentionally to compare data between bins.
extract_top_10_PFAMS <- function(table) {
# Choose top/bottom 10 most popular PFAMs. Create table with Phylum in wide format
  
  source_name <- deparse(substitute(table))

  hi_phylum <- table %>% 
    filter(label == 'highest') %>% 
    group_by(PFAM, Phylum) %>% 
    summarise(count = n(), .groups = 'drop') %>% 
    pivot_wider(names_from = Phylum, values_from = 'count', values_fill = list(count = 0))
  
  hi <- table %>%
    filter(label == 'highest') %>%
    group_by(across(all_of(c("PFAM", 'description')))) %>%  # Dynamic grouping
    distinct(PFAM, .keep_all = TRUE) %>%
    summarise(n_count = sum(n_count), mean_tAI = mean(tAI, na.rm = TRUE), .groups = 'drop') %>%
    arrange(desc(n_count)) %>%
    slice_head(n = 10) %>%
    mutate(label = "highest", source_df = source_name, additional = '')
  
  hi <- left_join(hi, hi_phylum, by = 'PFAM') %>% 
    replace_na(list(count = 0))

  
  
  lo_phylum <- table %>% 
    filter(label == 'lowest') %>% 
    group_by(PFAM, Phylum) %>% 
    summarise(count = n(), .groups = 'drop') %>% 
    pivot_wider(names_from = Phylum, values_from = 'count', values_fill = list(count = 0))

  lo <- table %>%
    filter(label == 'lowest') %>%
    group_by(across(all_of(c("PFAM", 'description')))) %>%  # Dynamic grouping
    distinct(PFAM, .keep_all = TRUE) %>%
    summarise(n_count = sum(n_count), mean_tAI = mean(tAI, na.rm = TRUE), .groups = 'drop') %>%
    arrange(desc(n_count)) %>%
    slice_head(n = 10) %>%
    mutate(label = "lowest", source_df = source_name, additional = '')

  
  lo <- left_join(lo, lo_phylum, by = 'PFAM') %>% 
    replace_na(list(count = 0))
  
  
  bind_rows(hi, lo)
}



count_enrichment <- function(original_df, extremum_df) {
  phylums <- c("Ascomycota", "Basidiomycota", "Entomophthoromycota", "Kickxellomycota", "Mortierellomycota", "Mucoromycota", "Chytridiomycota", "Microsporidia", "Monoblepharomycota", "Rozellomycota", "Glomeromycota", "Neocallimastigomycota")
  
  extremum_pfams <- unique(extremum_df$PFAM)
  original_df <- unique(original_df)
  extremum_df <- unique(extremum_df)
  
  counted_pfam_per_phylum <- original_df %>% 
    filter(PFAM %in% extremum_pfams) %>% 
    group_by(Phylum, PFAM) %>% 
    summarise(Count = n(), .groups='drop') %>% 
    group_by(Phylum) %>% 
    mutate(percentage = Count / sum(Count)) %>% 
    ungroup()
  
  #counts pfam occurences overall
  counted_pfam <- original_df %>%
    select(PFAM, prot_id) %>%
    filter(PFAM %in% extremum_pfams) %>%
    group_by(PFAM) %>%
    summarise(count = n()) %>%
    mutate(percent = count / sum(count)) %>%
    ungroup()
    


  extremum_long <- extremum_df %>%
    pivot_longer(
      cols = all_of(phylums),
      names_to = "Phylum",
      values_to = "Count" #there Count means Count_PFAMs_in_phylum
    ) %>%
    group_by(Phylum) %>%
    mutate(percentage = Count / sum(Count)) %>%
    ungroup()



  merged_df <- left_join(extremum_long, counted_pfam_per_phylum, by=c("PFAM", "Phylum")) %>%
    mutate(
      enrichment = ifelse(is.na(percentage.x) | percentage.y == 0, NA, percentage.x / percentage.y) #count.x are values in extremum df, count.y are values from original_df
    ) %>%
    rename(Count = Count.x) %>%
    select(PFAM, description, n_count, mean_tAI, Count,  label, source_df, additional, Phylum, enrichment) %>%
    pivot_wider(names_from = Phylum, values_from=c(Count, enrichment))

  result <- left_join(merged_df, counted_pfam, by = 'PFAM')
  result <- result%>%
    mutate(enrichment_overall = ((n_count / sum(n_count)) / percent))
  }


#```


#```{r create extreme tAI values PFAMS}
# tutaj tworze kilka tabel.
# ..._extreme - to te, które zawierają ileś % skrajnych białek. Tutaj mogą być zduplikowane białka, ponieważ niektóre białka zawierają kilka domen PFAM. Na wcześniejszym etapie rozdzielono je. Ułatwia to potem zliczanie liczby PFAMów co jest potrzebne w następnym etapie.
# hilo_... - Zlicza liczbę występień PFAMów we wszystkich białkach
# most_popular_PFAMS - wybiera z hilo... 10 najczęstszych (po 5 dla hi oraz lo) i łączy to dla wszystkich tabel hilo...


percent_to_get_extreme_tAI_values <- 10


#funkcja count_extreme_goterms  do poprawy
most_popular_PFAMS <- data.frame() #dataframe that store f.e. top 5 most common proteins with certain feature from df with top10% tAI proteins 

# Extract proteins with extreme tAI values (top and bottom 5%) from GOterms data
goterms_extreme <- get_extreme_tAI_values(GOterms_shorted, 'tAI', percent_to_get_extreme_tAI_values, 'percent') # Select top and bottom 5% tAI values
goterms_extreme <- goterms_extreme %>% 
  group_by(label, PFAM, GOterms, prot_id, tAI) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "") # Remove rows with missing PFAM annotations

# Count most frequent PFAM domains among extreme tAI proteins
hilo_goterms <- count_extreme_goterms(goterms_extreme, 10, FALSE)
hilo_goterms <- hilo_goterms %>% 
  group_by(label) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
write.table(hilo_goterms, file = paste0(main_path, "/hilo_goterms.csv"), sep = ",", row.names = FALSE, quote = FALSE)

# most_popular_PFAMS <- extract_top_10_PFAMS(hilo_goterms)


# Repeat the same procedure for proteins grouped by LCR length
LCR_length_extreme <- get_extreme_tAI_values(LCR_length, 'tAI', percent_to_get_extreme_tAI_values, 'percent')
LCR_length_extreme <- LCR_length_extreme %>% 
  group_by(label, PFAM, GOterms, prot_id, tAI) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")

#tabela ta daje nam informacje:
#W jakich białkach występują najczęstsze PFAM oraz z jakimi charakterystykami. W n_count podana została łączna liczba wystąpień PFAM w całym tym zbiorze. Duplikują się białka bo niektóre z nich mogą mieć kilka domen PFAM
hilo_LCR_len <- count_extreme_goterms(LCR_length_extreme, 10, FALSE, 'bin')
hilo_LCR_len <- hilo_LCR_len %>% 
  group_by(label) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
write.table(hilo_LCR_len, file = paste0(output_path, "hilo_LCR_len.csv"), sep = ",", row.names = FALSE, quote = FALSE)

most_popular_PFAMS <-extract_top_10_PFAMS(hilo_LCR_len)




# Analyze extreme tAI in relation to the number of LCRs
LCR_count_extreme <- get_extreme_tAI_values(LCR_number, 'tAI', percent_to_get_extreme_tAI_values, 'percent')
LCR_count_extreme <- LCR_count_extreme %>% 
  group_by(label, PFAM, GOterms, prot_id, tAI) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
hilo_LCR_presence <- count_extreme_goterms(LCR_count_extreme, 10, FALSE, 'bin')
hilo_LCR_presence <- hilo_LCR_presence %>% 
  group_by(label) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
write.table(hilo_LCR_presence, file = paste0(output_path, "LCR_count_extreme.csv"), sep = ",", row.names = FALSE, quote = FALSE)

most_popular_PFAMS <- bind_rows(most_popular_PFAMS, extract_top_10_PFAMS(hilo_LCR_presence))



# Analyze proteins with PFAMs associated with LCRs
# PFAM_LCR_extreme <- get_extreme_tAI_values(PFAM_LCR, 'tAI', 10, 'percent')
# PFAM_LCR_extreme <- PFAM_LCR_extreme %>% 
#   group_by(label, PFAM, GOterms, prot_id, tAI) %>%
#   filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
# hilo_PFAM_LCR <- count_extreme_goterms(PFAM_LCR_extreme, 10, FALSE)
# hilo_PFAM_LCR <- hilo_PFAM_LCR %>% 
#   group_by(label) %>%
#   filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
# write.table(hilo_PFAM_LCR, file = paste0(output_path, "hilo_PFAM_LCR.csv"), sep = ",", row.names = FALSE, quote = FALSE)
# 
# most_popular_PFAMS <- bind_rows(most_popular_PFAMS, extract_top_10_PFAMS(LCR_count_extreme))



# Evaluate impact of protein length on tAI extremes
prot_len_extreme <- get_extreme_tAI_values(protein_length, 'tAI', percent_to_get_extreme_tAI_values, 'percent')
prot_len_extreme <- prot_len_extreme %>% 
  group_by(label, PFAM, GOterms, prot_id, tAI) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
hilo_prot_len <- count_extreme_goterms(prot_len_extreme, 10, FALSE, 'bin')
hilo_prot_len <- hilo_prot_len %>%
  group_by(label) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
write.table(hilo_prot_len, file = paste0(output_path, "hilo_prot_len.csv"), sep = ",", row.names = FALSE, quote = FALSE)

most_popular_PFAMS <- bind_rows(most_popular_PFAMS, extract_top_10_PFAMS(hilo_prot_len))




# Analyze tAI extremes for proteins with or without signal peptides
signal_extreme <- get_extreme_tAI_values(signal_bolean, 'tAI', percent_to_get_extreme_tAI_values, 'percent')
signal_extreme <- signal_extreme %>% 
  group_by(label, PFAM, GOterms, prot_id, tAI) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
hilo_signal <- count_extreme_goterms(signal_extreme, 10, FALSE)
hilo_signal <- hilo_signal %>% 
  group_by(label) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
write.table(hilo_signal, file = paste0(output_path, "hilo_signal.csv"), sep = ",", row.names = FALSE, quote = FALSE)

most_popular_PFAMS <- bind_rows(most_popular_PFAMS, extract_top_10_PFAMS(hilo_signal))



# Examine proteins with/without transmembrane regions (presence)
TM_pres_extreme <- get_extreme_tAI_values(TM_count, 'tAI', percent_to_get_extreme_tAI_values, 'percent')
TM_pres_extreme <- TM_pres_extreme %>% 
  group_by(label, PFAM, GOterms, prot_id, tAI) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
hilo_TM_presence <- count_extreme_goterms(TM_pres_extreme, 10, FALSE, 'bin')
hilo_TM_presence <- hilo_TM_presence %>% 
  group_by(label) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
write.table(hilo_TM_presence, file = paste0(output_path, "hilo_TM_presence.csv"), sep = ",", row.names = FALSE, quote = FALSE)

most_popular_PFAMS <- bind_rows(most_popular_PFAMS, extract_top_10_PFAMS(hilo_TM_presence))




# Count of transmembrane domains vs tAI extremes
TM_count_extreme <- get_extreme_tAI_values(TM_count, 'tAI', percent_to_get_extreme_tAI_values, 'percent')
TM_count_extreme <- TM_count_extreme %>% 
  group_by(label, PFAM, GOterms, prot_id, tAI) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
hilo_TM_count <- count_extreme_goterms(TM_count_extreme, 10, FALSE, 'bin')
hilo_TM_count <- hilo_TM_count %>% 
  group_by(label) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
write.table(hilo_TM_count, file = paste0(output_path, "hilo_TM_count.csv"), sep = ",", row.names = FALSE, quote = FALSE)

most_popular_PFAMS <- bind_rows(most_popular_PFAMS, extract_top_10_PFAMS(hilo_TM_count))




# Analyze all variables together for extreme tAI values
all_values_extreme <- get_extreme_tAI_values(all_values_table, 'tAI', percent_to_get_extreme_tAI_values, 'percent')
all_values_extreme <- all_values_extreme %>% 
  group_by(label, PFAM, GOterms, prot_id, acc_id, tAI) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
hilo_all_values_count <- count_extreme_goterms(all_values_extreme, 10, FALSE)
hilo_all_values_count <- hilo_all_values_count %>% 
  group_by(label) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
write.table(hilo_all_values_count, file = paste0(output_path, "hilo_all_values_count.csv"), sep = ",", row.names = FALSE, quote = FALSE)


# Focus only on proteins with PFAMs in the all-values dataset
all_values_only_PFAM_extreme <- get_extreme_tAI_values(all_values_only_PFAMs, 'tAI', percent_to_get_extreme_tAI_values, 'percent')
all_values_only_PFAM_extreme <- all_values_only_PFAM_extreme %>% 
  group_by(label, PFAM, GOterms, prot_id, acc_id, tAI) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
hilo_all_values_ONLY_PFAMs_count <- count_extreme_goterms(all_values_only_PFAM_extreme, 10, FALSE)
hilo_all_values_ONLY_PFAMs_count <- hilo_all_values_ONLY_PFAMs_count %>% 
  group_by(label) %>%
  filter(!is.na(PFAM), !is.nan(PFAM), PFAM != "")
write.table(hilo_all_values_ONLY_PFAMs_count, file = paste0(output_path, "hilo_all_values_ONLY_PFAMs_count.csv"), sep = ",", row.names = FALSE, quote = FALSE)

most_popular_PFAMS <- bind_rows(most_popular_PFAMS, extract_top_10_PFAMS(hilo_all_values_ONLY_PFAMs_count))
most_popular_PFAMS <- most_popular_PFAMS %>% 
  mutate_all(~ replace(., is.na(.), 0))
most_popular_PFAMS <- count_enrichment(all_values_only_PFAMs, most_popular_PFAMS)
  

# Merge all previously created datasets
hilo_merged <- rbind(hilo_goterms, hilo_LCR_len, hilo_LCR_presence, hilo_prot_len, hilo_signal, hilo_TM_count, hilo_TM_presence)
hilo_all_values_merged <- rbind(all_values_only_PFAM_extreme, all_values_extreme)



# Save merged data
write.table(hilo_merged, file = paste0(output_path, "hilo_merged.csv"), sep = ",", row.names = FALSE, quote = FALSE)
write.table(hilo_all_values_merged, file = paste0(output_path, "hilo_all_values_merged.csv"), sep = ",", row.names = FALSE, quote = FALSE)

write.table(most_popular_PFAMS, file = paste0(output_path, "most_popular_PFAMs_in_tAI_tails.csv"), sep = "," , row.names = FALSE, quote = FALSE)

write.csv(most_popular_PFAMS, file = "most_popular_PFAMs_in_tAI_tails.csv", row.names = FALSE)


#```

# #```{r functions to secnd part}
# 
# #Calculates summary statistics (either percentage of valid entries or mean values) for selected columns in a data frame.
# count_parameters <- function(df, columns, variable) {
#     data <- df[, columns]  # Select relevant columns
#     
#     if (variable == "percent") {
#         total <- nrow(data)  # Total row count
#         count <- colSums(!is.na(data) & data != "" & data != 0)  # Count valid entries
#         percent <- (count / total) * 100  # Calculate percentage
#         
#         return(list(total = total, count = count, percent = percent))  # Ensure function returns the result
#     } else if (variable == 'mean') {
#         mean_vals <- sapply(data, function(x) mean(x[!is.na(x) & x != 0], na.rm = TRUE))
# 
#         return(mean_vals)
#     }
# }
# 
# 
# 
# 
# #Combines a list of results (from count_parameters) into a single, long-format data frame for analysis or visualization.
# merge_results <- function(dfs_list, variable) {
#   result <- do.call(rbind, lapply(names(dfs_list), function(name) {
#     res <- dfs_list[[name]]
#     
#     if (variable == 'percent') {
#       data.frame(
#         dataset = name,
#         variable = names(res$percent),
#         percent = as.numeric(res$percent),
#         count = as.numeric(res$count),
#         total = res$total
#       )
#     } else if (variable == 'mean') {
#       data.frame(
#         dataset = name,
#         variable = names(res),
#         mean = as.numeric(res)
#       )
#     }
#   }))
#   
#   return(result)
# }
# 
# #```
# 
# #```{r ten thousand top/bottom proteins exploration}
# 
# tenk_hi<- all_values_table %>% 
#   arrange(tAI) %>% 
#   head(10000)
# tenk_lo  <- all_values_table %>% 
#   arrange(tAI) %>% 
#   tail(10000)
# 
# tmp <- count_parameters(tenk_hi, c('domain', 'signal'), 'percent')
# 
# 
# 
# #```


