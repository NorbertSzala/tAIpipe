#```{r TO DO}
####### DONE #########
# split big_df to 1000hi and 1000lo datasets. count how many is there proteins with: a) LCR, b) what is LCR length, c) prot len, d) signal peptides etc.
# Compare it between has pfams, doesnt have pfams
# Compare it to main dataset



####### TO DO ########
#make 20000 long datafram (half for highest, half for lowest tAI)
# Check what type of proteins make up dataframe/ - if they consist LCR, TM etc...


#```



#```{r libraries}
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggExtra)
library(ggthemes)

#```



#```{r functions}

get_extreme_tAI_values <- function(table, column, extremum, variable) {
  #' Get Extreme tAI Values
  #' This function returns a dataframe with proteins that have extreme tAI values from top bins 
  #' (the lowest and highest) depending on the bins. The 'extremum' parameter 
  #' is an integer representing the top 'extremum' percentage of values.
  #'
  #' @param table A data frame containing the data.
  #' @param column The name of the column containing tAI values.
  #' @param extremum The top 'extremum' percentage of values or number of bins
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




#Calculates summary statistics (either percentage of valid entries or mean values) for selected columns in a data frame.
count_parameters <- function(df, columns, variable) {
    data <- df[, columns]  # Select relevant columns
    
    if (variable == "percent") {
        total <- nrow(data)  # Total row count
        count <- colSums(!is.na(data) & data != "" & data != 0)  # Count valid entries
        percent <- (count / total) * 100  # Calculate percentage

        
        return(list(total = total, count = count, percent = percent))  # Ensure function returns the result
    } else if (variable == 'mean') {
        mean_vals <- sapply(data, function(x) mean(x[!is.na(x) & x != 0], na.rm = TRUE))

        return(mean_vals)
    }
}




#Combines a list of results (from count_parameters) into a single, long-format data frame for analysis or visualization.
merge_results <- function(dfs_list, variable) {
  result <- do.call(rbind, lapply(names(dfs_list), function(name) {
    res <- dfs_list[[name]]
    
    if (variable == 'percent') {
      data.frame(
        dataset = name,
        variable = names(res$percent),
        percent = as.numeric(res$percent),
        count = as.numeric(res$count),
        total = res$total
      )
    } else if (variable == 'mean') {
      data.frame(
        dataset = name,
        variable = names(res),
        mean = as.numeric(res)
      )
    }
  }))
  
  return(result)
}


#```


#```{r read data}
# local path
main_path <- '/home/norbert/IBB/skrypty/Main_Dataset_splitted/data/'
output_path<- '/home/norbert/IBB/skrypty/Main_Dataset_splitted/plots_top_and_bottom_tAI/'
plots_path <- '/home/norbert/IBB/skrypty/Main_Dataset_splitted/plots_all_vs_top_bottom/'


#server mycelia path
# main_path <- '/home/norbert_s/tAI/all_chunks/Main_Dataset_splitted/data/'
# output_path<-'/home/norbert_s/tAI/all_chunks/Main_Dataset_splitted/plots_extensive



PFAM_acc_de <- read.csv(paste0(main_path, 'PFAM_acc_de.csv'), header = FALSE, col.names = c('accession', 'description'))

#its mostly the sam what all_values_table known from previous R scripts
big_df <- read.csv(paste0(main_path, 'shorted_main_dataset.tsv'), sep = '\t', header = FALSE, col.names = c('acc_id', 'prot_id','prot_len', 'domain_presence', 'TM_presence', 'TM_length', 'signal_presence', 'LCR_presence', 'LCR_length', 'PFAM_LCR', 'PFAM', 'GOterms', 'tAI'))

#split PFAM column to make each PFAM id in unique row. The same for goterms
# big_df <- all_values_table %>% 
  # separate_rows(PFAM, sep = " ") %>%
  # separate_rows(GOterms, sep = '\\|')
#move version symbol from PFAM ID
big_df$PFAM <- sub("\\..*$", "", big_df$PFAM)
big_df$PFAM_presence <- ifelse(big_df$PFAM != '', 1, 0)
big_df <- left_join(big_df, PFAM_acc_de, by = c('PFAM' = 'accession'))
#all values without labeling based on tAI

#top 10% proteins with tAI
med_df <- get_extreme_tAI_values(big_df, 'tAI', 10, 'percent')
med_hi <- med_df %>% 
    filter(label == 'highest')
med_lo <- med_df %>% 
    filter(label == 'lowest')

#top 1% proteins with tAI
small_df <- get_extreme_tAI_values(big_df, 'tAI', 1, 'percent')
small_hi <- med_df %>% 
    filter(label == 'highest')
small_lo <- med_df %>% 
    filter(label == 'lowest')


#```




#```{r data frames}

dfs <- list(big_df, med_hi, med_lo, small_hi, small_lo)
names(dfs) <- c("big_df", "med_hi", "med_lo", "small_hi", "small_lo")
cols_categorical <- c("domain_presence", "TM_presence", "signal_presence", "LCR_presence")
cols_continuous <-  c("prot_len", "TM_length", "LCR_length")




categ_res <- lapply(dfs, function(df) {count_parameters(df, cols_categorical, 'percent')})
cont_res <- lapply(dfs, function(df) {count_parameters(df, cols_continuous, 'mean')})





categ_df <- merge_results(categ_res, 'percent')
cont_df <- merge_results(cont_res, 'mean')
#```



#```{r plots}

#Wykres pokazuje jaki procent białek w danym zbiorze posiada cechę (np. obecność peptydów sygnałowych). 
categ_plot <- ggplot(categ_df, aes(x = variable, y = percent, fill = dataset)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  theme_stata() +
  theme(
    axis.text.x = element_text(face = 'bold', hjust = 1),
    axis.text.y = element_text(face = 'bold')
  ) +
    geom_text(
        aes(label = count),
        position = position_dodge(width = 0.9),
        vjust = -0.5,
        size = 3,
        color = 'black'
    )+
  labs(
    x = '',
    y = "Percent",
    title = "Comparing the Percentage of Feature Presence Across Different Datasets",
    fill = "Dataset"
  ) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_fill_brewer(palette = "Set1", labels = c('All proteins', 'Top 10%', 'Bottom 10%', 
               'Top 1%', 'Bottom 10%'))
categ_plot

ggsave(paste0(plots_path, 'categorical_variables_vs_all_proteins.png'), plot=categ_plot, width = 10, height = 8, dpi=300)




#count differences between whole proteins set and others
diff_means <- cont_df %>% 
    group_by(variable) %>% 
    mutate(mean = round(mean, 2),
           diff = round(mean - mean[dataset == 'big_df'], 2)
           ) %>% 
    ungroup() %>% 
    pull(diff)

#Wykres pokazuje jaki procent białek w danym zbiorze posiada cechę (np. obecność peptydów sygnałowych). 
cont_plot <- ggplot(cont_df, aes(x = variable, y = mean, fill = dataset)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  theme_stata() +
  theme(
    axis.text.x = element_text(face = 'bold', hjust = 1),
    axis.text.y = element_text(face = 'bold')
  ) +
    geom_text(
        aes(label = diff_means),
        position = position_dodge(width = 0.9),
        vjust = -0.5,
        size = 3,
        color = 'black'
    )+
  labs(
    x = '',
    y = "Mean length (AA)",
    title = "Comparing the mean lengths of Features Across Different Datasets",
    fill = "Dataset"
  ) +
    scale_fill_brewer(palette = "Set1", labels = c('All proteins', 'Top 10%', 'Bottom 10%', 
               'Top 1%', 'Bottom 10%'))
cont_plot

ggsave(paste0(plots_path, 'continuous_variables_vs_all_proteins.png'), plot=cont_plot, width = 10, height = 8, dpi=300)



# how to read theese plots?
#Protein domains are found in ~65% of all proteins. In contrast, among proteins with tAI values in the top 10%, domains are present in ~75% of the cases.
##```
#```

