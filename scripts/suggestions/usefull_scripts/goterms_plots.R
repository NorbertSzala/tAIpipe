#```{r TODO}
# Przekodować sekwencje LCR (aa) na nukleotydy i porownac  wartosci tego tAI z bez LCR
# Czytać
#```



#```{r libraries}
library(ggplot2)
library(ggthemes)
library(viridis)
library(dplyr)
library(RColorBrewer)
library(grDevices)
library(dunn.test)
library(scales)
library(tidyr)
library(ggExtra)
library(nortest)
library(cowplot)
library(stringr)
#```

#```{r paths}
main_path <- '/home/norbert/IBB/skrypty/Main_Dataset_splitted/data/'

chosen_GOterms <- read.csv(paste0(main_path, 'chosen_GOterms.tsv'), sep='\t', header=FALSE, col.names = c('ass_id', 'prot_id','prot_len', 'domain', 'TM_count', 'TM_length', 'signal', 'LCR_number', 'LCR_length', 'PFAM_LCR', 'PFAM', 'GOterms', 'tAI'))

output_path<-'/home/norbert/IBB/skrypty/Main_Dataset_splitted/plots_goterms/'
#```

#```{r read data - mycelia paths}
main_path <- '/home/norbert_s/tAI/all_chunks/Main_Dataset_splitted/data/'

chosen_GOterms <- read.csv(paste0(main_path, 'chosen_GOterms.tsv'), sep='\t', header=FALSE, col.names = c('ass_id', 'prot_id','prot_len', 'domain', 'TM_count', 'TM_length', 'signal', 'LCR_number', 'LCR_length', 'PFAM_LCR', 'PFAM', 'GOterms', 'tAI'))

output_path<-'/home/norbert_s/tAI/all_chunks/Main_Dataset_splitted/goterms_plots/'
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
    ad_test <- ad.test(var1) #Anderson Darling test to check normality distribution
    
    if (ad_test$p.value < 0.05) {
      # ANOVA test for parametric values with normal distribution
      test = 'ANOVA p-value: '
      anova_results <- summary(aov(var1 ~ var2, data = data))
      p_value <- anova_results[[1]]$'Pr(>F)'[1]
    } else {
      # Kruskal-Wallis test for parametric values and lack of normal distribution
      test = 'Kruskal-Wallis test p-value: '
      kw_test <- kruskal.test(var1 ~ var2, data = data)
      p_value <- kw_test$p.value
    }
  } else {
    ad_test <- ad.test(var1)
    
    if (ad_test$p.value < 0.05) {
      # Mann-Whitney test #non parametrical test with normal distribution
      test = 'Mann-Whitney test p-value: '
      mann_whitney_test <- wilcox.test(var1 ~ var2, data = data)
      p_value <- mann_whitney_test$p.value
    } else {
      # T-test for not normal test not parametrical
      test = 'T-student test p-value: '
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
#```


#```{r split data and make long table}

domain_GOterms <- chosen_GOterms %>% 
  mutate(
    RNA = case_when(
      str_detect(GOterms, 'RNA') & domain == 1 ~ 'RNA_DOMP',
      TRUE ~ 'RNA_DOMA'),
    DNA = case_when(
      str_detect(GOterms, 'DNA') & domain == 1 ~ 'DNA_DOMP',
      TRUE ~ 'DNA_DOMA'),
    Mito = case_when(
      str_detect(GOterms, 'mitochondrial') & domain == 1 ~ 'Mito_DOMP',
      TRUE ~ 'Mito_DOMA'),
    Prot = case_when(
      str_detect(GOterms, 'protein') & domain == 1 ~ 'Prot_DOMP',
      TRUE ~ 'Prot_DOMA'),
    Lipid = case_when(
      str_detect(GOterms, 'lipid') & domain == 1 ~ 'Lip_DOMP',
      TRUE ~ 'Lip_DOMA'),
    ) %>% 
  select('prot_id', 'tAI', 'RNA', 'DNA', 'Mito', 'Prot', 'Lipid') %>% 
  pivot_longer(cols=c('RNA', 'DNA', 'Mito', 'Prot', 'Lipid'),
               names_to = 'GO_term',
               values_to = 'Presence')

TM_GOterms <- chosen_GOterms %>% 
  mutate(
    RNA = case_when(
      str_detect(GOterms, 'RNA') & TM_count == 1 ~ 'RNA_TM1',
      str_detect(GOterms, 'RNA') & TM_count > 1 ~ 'RNA_TM>1',
      TRUE ~ 'RNA_TMA'),
    DNA = case_when(
      str_detect(GOterms, 'DNA') & TM_count == 1 ~ 'DNA_TM1',
      str_detect(GOterms, 'DNA') & TM_count > 1 ~ 'DNA_TM>1',
      TRUE ~ 'DNA_TMA'),
    Mito = case_when(
      str_detect(GOterms, 'mitochondrial') & TM_count == 1 ~ 'Mito_TM1',
      str_detect(GOterms, 'Mitochondrial') & TM_count > 1 ~ 'Mito_TM>1',
      TRUE ~ 'Mito_TMA'),
    Prot = case_when(
      str_detect(GOterms, 'protein') & TM_count == 1 ~ 'Prot_TM1',
      str_detect(GOterms, 'protein') & TM_count > 1 ~ 'Prot_TM>1',
      TRUE ~ 'Prot_TMA'),
    Lipid = case_when(
      str_detect(GOterms, 'lipid') & TM_count == 1 ~ 'Lip_TM1',
      str_detect(GOterms, 'lipid') & TM_count > 1 ~ 'Lip_TM>1',
      TRUE ~ 'Lip_TMA')
  ) %>% 
  select(prot_id, tAI, RNA, DNA, Mito, Prot, Lipid) %>% 
  pivot_longer(cols = c(RNA, DNA, Mito, Prot, Lipid),
               names_to = 'GO_term',
               values_to = 'Presence')



chosen_GOterms$TM_length <- as.numeric(chosen_GOterms$TM_length)
max_len <- max(chosen_GOterms$TM_length, na.rm = TRUE)
chosen_GOterms$TM_Len_bin <- cut(
  chosen_GOterms$TM_length,
  breaks = c(-Inf, 1, 0.25 * max_len, 0.5 * max_len, 0.75 * max_len, max_len),
  labels = c("0", "Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE,
  right = FALSE 
)

TM_len_GOterms <- chosen_GOterms %>%
  mutate(
    RNA = case_when(
    str_detect(GOterms, 'RNA') & domain == 1 & TM_Len_bin == '0' ~ 'RNA_TMP_0',
    str_detect(GOterms, 'RNA') & domain == 1 & TM_Len_bin == 'Q1' ~ 'RNA_TMP_Q1',
    str_detect(GOterms, 'RNA') & domain == 1 & TM_Len_bin == 'Q2' ~ 'RNA_TMP_Q2',
    str_detect(GOterms, 'RNA') & domain == 1 & TM_Len_bin == 'Q3' ~ 'RNA_TMP_Q3',
    str_detect(GOterms, 'RNA') & domain == 1 & TM_Len_bin == 'Q4' ~ 'RNA_TMP_Q4',
    str_detect(GOterms, 'RNA') & domain == 0 & TM_Len_bin == '0' ~ 'RNA_TMA_0',
    str_detect(GOterms, 'RNA') & domain == 0 & TM_Len_bin == 'Q1' ~ 'RNA_TMA_Q1',
    str_detect(GOterms, 'RNA') & domain == 0 & TM_Len_bin == 'Q2' ~ 'RNA_TMA_Q2',
    str_detect(GOterms, 'RNA') & domain == 0 & TM_Len_bin == 'Q3' ~ 'RNA_TMA_Q3',
    str_detect(GOterms, 'RNA') & domain == 0 & TM_Len_bin == 'Q4' ~ 'RNA_TMA_Q4'),
    
    DNA = case_when(
    str_detect(GOterms, 'DNA') & domain == 1 & TM_Len_bin == '0' ~ 'DNA_TMP_0',
    str_detect(GOterms, 'DNA') & domain == 1 & TM_Len_bin == 'Q1' ~ 'DNA_TMP_Q1',
    str_detect(GOterms, 'DNA') & domain == 1 & TM_Len_bin == 'Q2' ~ 'DNA_TMP_Q2',
    str_detect(GOterms, 'DNA') & domain == 1 & TM_Len_bin == 'Q3' ~ 'DNA_TMP_Q3',
    str_detect(GOterms, 'DNA') & domain == 1 & TM_Len_bin == 'Q4' ~ 'DNA_TMP_Q4',
    str_detect(GOterms, 'DNA') & domain == 0 & TM_Len_bin == '0' ~ 'DNA_TMA_0',
    str_detect(GOterms, 'DNA') & domain == 0 & TM_Len_bin == 'Q1' ~ 'DNA_TMA_Q1',
    str_detect(GOterms, 'DNA') & domain == 0 & TM_Len_bin == 'Q2' ~ 'DNA_TMA_Q2',
    str_detect(GOterms, 'DNA') & domain == 0 & TM_Len_bin == 'Q3' ~ 'DNA_TMA_Q3',
    str_detect(GOterms, 'DNA') & domain == 0 & TM_Len_bin == 'Q4' ~ 'DNA_TMA_Q4'),
    
    Mito = case_when(
    str_detect(GOterms, 'mitochondrial') & domain == 1 & TM_Len_bin == '0' ~ 'Mito_TMP_0',
    str_detect(GOterms, 'mitochondrial') & domain == 1 & TM_Len_bin == 'Q1' ~ 'Mito_TMP_Q1',
    str_detect(GOterms, 'mitochondrial') & domain == 1 & TM_Len_bin == 'Q2' ~ 'Mito_TMP_Q2',
    str_detect(GOterms, 'mitochondrial') & domain == 1 & TM_Len_bin == 'Q3' ~ 'Mito_TMP_Q3',
    str_detect(GOterms, 'mitochondrial') & domain == 1 & TM_Len_bin == 'Q4' ~ 'Mito_TMP_Q4',
    str_detect(GOterms, 'mitochondrial') & domain == 0 & TM_Len_bin == '0' ~ 'Mito_TMA_0',
    str_detect(GOterms, 'mitochondrial') & domain == 0 & TM_Len_bin == 'Q1' ~ 'Mito_TMA_Q1',
    str_detect(GOterms, 'mitochondrial') & domain == 0 & TM_Len_bin == 'Q2' ~ 'Mito_TMA_Q2',
    str_detect(GOterms, 'mitochondrial') & domain == 0 & TM_Len_bin == 'Q3' ~ 'Mito_TMA_Q3',
    str_detect(GOterms, 'mitochondrial') & domain == 0 & TM_Len_bin == 'Q4' ~ 'Mito_TMA_Q4'),
    
    Prot = case_when(
    str_detect(GOterms, 'protein') & domain == 1 & TM_Len_bin == '0' ~ 'Prot_TMP_0',
    str_detect(GOterms, 'protein') & domain == 1 & TM_Len_bin == 'Q1' ~ 'Prot_TMP_Q1',
    str_detect(GOterms, 'protein') & domain == 1 & TM_Len_bin == 'Q2' ~ 'Prot_TMP_Q2',
    str_detect(GOterms, 'protein') & domain == 1 & TM_Len_bin == 'Q3' ~ 'Prot_TMP_Q3',
    str_detect(GOterms, 'protein') & domain == 1 & TM_Len_bin == 'Q4' ~ 'Prot_TMP_Q4',
    str_detect(GOterms, 'protein') & domain == 0 & TM_Len_bin == '0' ~ 'Prot_TMA_0',
    str_detect(GOterms, 'protein') & domain == 0 & TM_Len_bin == 'Q1' ~ 'Prot_TMA_Q1',
    str_detect(GOterms, 'protein') & domain == 0 & TM_Len_bin == 'Q2' ~ 'Prot_TMA_Q2',
    str_detect(GOterms, 'protein') & domain == 0 & TM_Len_bin == 'Q3' ~ 'Prot_TMA_Q3',
    str_detect(GOterms, 'protein') & domain == 0 & TM_Len_bin == 'Q4' ~ 'Prot_TMA_Q4'),
    
    Lipid = case_when(
    str_detect(GOterms, 'lipid') & domain == 1 & TM_Len_bin == '0' ~ 'Lip_TMP_0',
    str_detect(GOterms, 'lipid') & domain == 1 & TM_Len_bin == 'Q1' ~ 'Lip_TMP_Q1',
    str_detect(GOterms, 'lipid') & domain == 1 & TM_Len_bin == 'Q2' ~ 'Lip_TMP_Q2',
    str_detect(GOterms, 'lipid') & domain == 1 & TM_Len_bin == 'Q3' ~ 'Lip_TMP_Q3',
    str_detect(GOterms, 'lipid') & domain == 1 & TM_Len_bin == 'Q4' ~ 'Lip_TMP_Q4',
    str_detect(GOterms, 'lipid') & domain == 0 & TM_Len_bin == '0' ~ 'Lip_TMA_0',
    str_detect(GOterms, 'lipid') & domain == 0 & TM_Len_bin == 'Q1' ~ 'Lip_TMA_Q1',
    str_detect(GOterms, 'lipid') & domain == 0 & TM_Len_bin == 'Q2' ~ 'Lip_TMA_Q2',
    str_detect(GOterms, 'lipid') & domain == 0 & TM_Len_bin == 'Q3' ~ 'Lip_TMA_Q3',
    str_detect(GOterms, 'lipid') & domain == 0 & TM_Len_bin == 'Q4' ~ 'Lip_TMA_Q4')) %>% 
  select('prot_id', 'tAI', 'RNA', 'DNA', 'Mito', 'Prot', 'Lipid') %>%  # removed extra comma here
  pivot_longer(cols=c('RNA', 'DNA', 'Mito', 'Prot', 'Lipid'),
               names_to = 'GO_term',
               values_to = 'Presence') %>% 
  filter(!is.na(Presence))



signal_GOterms <- chosen_GOterms %>% 
  mutate(
    RNA = case_when(
      str_detect(GOterms, 'RNA') & signal == 1 ~ 'RNA_SIGP',
      TRUE ~ 'RNA_SIGA'),
    DNA = case_when(
      str_detect(GOterms, 'DNA') & signal == 1 ~ 'DNA_SIGP',
      TRUE ~ 'DNA_SIGA'),
    Mito = case_when(
      str_detect(GOterms, 'mitochondrial') & signal == 1 ~ 'Mito_SIGP',
      TRUE ~ 'Mito_SIGA'),
    Prot = case_when(
      str_detect(GOterms, 'protein') & signal == 1 ~ 'Prot_SIGP',
      TRUE ~ 'Prot_SIGA'),
    Lipid = case_when(
      str_detect(GOterms, 'lipid') & signal == 1 ~ 'Lip_SIGP',
      TRUE ~ 'Lip_SIGA'),
    ) %>% 
  select('prot_id', 'tAI', 'RNA', 'DNA', 'Mito', 'Prot', 'Lipid') %>% 
  pivot_longer(cols=c('RNA', 'DNA', 'Mito', 'Prot', 'Lipid'),
               names_to = 'GO_term',
               values_to = 'Presence')

LCR_GOterms <- chosen_GOterms %>% 
  mutate(
    RNA = case_when(
      str_detect(GOterms, 'RNA') & LCR_number == 1 ~ 'RNA_LCRP',
      str_detect(GOterms, 'RNA') & LCR_number > 1 ~ 'RNA_TM>1',
      TRUE ~ 'RNA_LCRA'),
    DNA = case_when(
      str_detect(GOterms, 'DNA') & LCR_number == 1 ~ 'DNA_LCRP',
      str_detect(GOterms, 'DNA') & LCR_number > 1 ~ 'DNA_LCRP>1',
      TRUE ~ 'DNA_LCRA'),
    Mito = case_when(
      str_detect(GOterms, 'mitochondrial') & LCR_number == 1 ~ 'Mito_LCRP',
      str_detect(GOterms, 'mitochondrial') & LCR_number > 1 ~ 'Mito_LCRP>1',
      TRUE ~ 'Mito_LCRA'),
    Prot = case_when(
      str_detect(GOterms, 'protein') & LCR_number == 1 ~ 'Prot_LCRP',
      str_detect(GOterms, 'protein') & LCR_number > 1 ~ 'Prot_LCRP>1',
      TRUE ~ 'Prot_LCRA'),
    Lipid = case_when(
      str_detect(GOterms, 'lipid') & LCR_number == 1 ~ 'Lip_LCRP',
      str_detect(GOterms, 'lipid') & LCR_number > 1 ~ 'Lip_LCRP>1',
      TRUE ~ 'Lip_LCRA'),
    ) %>% 
  select('prot_id', 'tAI', 'RNA', 'DNA', 'Mito', 'Prot', 'Lipid') %>% 
  pivot_longer(cols=c('RNA', 'DNA', 'Mito', 'Prot', 'Lipid'),
               names_to = 'GO_term',
               values_to = 'Presence')



chosen_GOterms$LCR_Len_bin <- as.numeric(chosen_GOterms$LCR_length)
max_len <- max(chosen_GOterms$LCR_length, na.rm = TRUE)
chosen_GOterms$LCR_Len_bin <- cut(
  chosen_GOterms$LCR_length,
  breaks = c(-Inf, 1, 0.25 * max_len, 0.5 * max_len, 0.75 * max_len, max_len),
  labels = c("0", "Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE,
  right = FALSE 
)




chosen_GOterms$LCR_length <- as.numeric(chosen_GOterms$LCR_length)
max_len <- max(chosen_GOterms$LCR_length, na.rm = TRUE)
chosen_GOterms$LCR_Len_bin <- cut(
  chosen_GOterms$LCR_length,
  breaks = c(-Inf, 1, 0.25 * max_len, 0.5 * max_len, 0.75 * max_len, max_len),
  labels = c("0", "Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE,
  right = FALSE 
)

LCR_len_GOterms <- chosen_GOterms %>%
  mutate(
    LCR_length = as.numeric(LCR_length),
    RNA = case_when(
    str_detect(GOterms, 'RNA') & domain == 1 & LCR_Len_bin == '0' ~ 'RNA_LCRP_0',
    str_detect(GOterms, 'RNA') & domain == 1 & LCR_Len_bin == 'Q1' ~ 'RNA_LCRP_Q1',
    str_detect(GOterms, 'RNA') & domain == 1 & LCR_Len_bin == 'Q2' ~ 'RNA_LCRP_Q2',
    str_detect(GOterms, 'RNA') & domain == 1 & LCR_Len_bin == 'Q3' ~ 'RNA_LCRP_Q3',
    str_detect(GOterms, 'RNA') & domain == 1 & LCR_Len_bin == 'Q4' ~ 'RNA_LCRP_Q4',
    str_detect(GOterms, 'RNA') & domain == 0 & LCR_Len_bin == '0' ~ 'RNA_LCRA_0',
    str_detect(GOterms, 'RNA') & domain == 0 & LCR_Len_bin == 'Q1' ~ 'RNA_LCRA_Q1',
    str_detect(GOterms, 'RNA') & domain == 0 & LCR_Len_bin == 'Q2' ~ 'RNA_LCRA_Q2',
    str_detect(GOterms, 'RNA') & domain == 0 & LCR_Len_bin == 'Q3' ~ 'RNA_LCRA_Q3',
    str_detect(GOterms, 'RNA') & domain == 0 & LCR_Len_bin == 'Q4' ~ 'RNA_LCRA_Q4'),
    
    DNA = case_when(
    str_detect(GOterms, 'DNA') & domain == 1 & LCR_Len_bin == '0' ~ 'DNA_LCRP_0',
    str_detect(GOterms, 'DNA') & domain == 1 & LCR_Len_bin == 'Q1' ~ 'DNA_LCRP_Q1',
    str_detect(GOterms, 'DNA') & domain == 1 & LCR_Len_bin == 'Q2' ~ 'DNA_LCRP_Q2',
    str_detect(GOterms, 'DNA') & domain == 1 & LCR_Len_bin == 'Q3' ~ 'DNA_LCRP_Q3',
    str_detect(GOterms, 'DNA') & domain == 1 & LCR_Len_bin == 'Q4' ~ 'DNA_LCRP_Q4',
    str_detect(GOterms, 'DNA') & domain == 0 & LCR_Len_bin == '0' ~ 'DNA_LCRA_0',
    str_detect(GOterms, 'DNA') & domain == 0 & LCR_Len_bin == 'Q1' ~ 'DNA_LCRA_Q1',
    str_detect(GOterms, 'DNA') & domain == 0 & LCR_Len_bin == 'Q2' ~ 'DNA_LCRA_Q2',
    str_detect(GOterms, 'DNA') & domain == 0 & LCR_Len_bin == 'Q3' ~ 'DNA_LCRA_Q3',
    str_detect(GOterms, 'DNA') & domain == 0 & LCR_Len_bin == 'Q4' ~ 'DNA_LCRA_Q4'),
    
    Mito = case_when(
    str_detect(GOterms, 'mitochondrial') & domain == 1 & LCR_Len_bin == '0' ~ 'Mito_LCRP_0',
    str_detect(GOterms, 'mitochondrial') & domain == 1 & LCR_Len_bin == 'Q1' ~ 'Mito_LCRP_Q1',
    str_detect(GOterms, 'mitochondrial') & domain == 1 & LCR_Len_bin == 'Q2' ~ 'Mito_LCRP_Q2',
    str_detect(GOterms, 'mitochondrial') & domain == 1 & LCR_Len_bin == 'Q3' ~ 'Mito_LCRP_Q3',
    str_detect(GOterms, 'mitochondrial') & domain == 1 & LCR_Len_bin == 'Q4' ~ 'Mito_LCRP_Q4',
    str_detect(GOterms, 'mitochondrial') & domain == 0 & LCR_Len_bin == '0' ~ 'Mito_LCRA_0',
    str_detect(GOterms, 'mitochondrial') & domain == 0 & LCR_Len_bin == 'Q1' ~ 'Mito_LCRA_Q1',
    str_detect(GOterms, 'mitochondrial') & domain == 0 & LCR_Len_bin == 'Q2' ~ 'Mito_LCRA_Q2',
    str_detect(GOterms, 'mitochondrial') & domain == 0 & LCR_Len_bin == 'Q3' ~ 'Mito_LCRA_Q3',
    str_detect(GOterms, 'mitochondrial') & domain == 0 & LCR_Len_bin == 'Q4' ~ 'Mito_LCRA_Q4'),
    
    Prot = case_when(
    str_detect(GOterms, 'protein') & domain == 1 & LCR_Len_bin == '0' ~ 'Prot_LCRP_0',
    str_detect(GOterms, 'protein') & domain == 1 & LCR_Len_bin == 'Q1' ~ 'Prot_LCRP_Q1',
    str_detect(GOterms, 'protein') & domain == 1 & LCR_Len_bin == 'Q2' ~ 'Prot_LCRP_Q2',
    str_detect(GOterms, 'protein') & domain == 1 & LCR_Len_bin == 'Q3' ~ 'Prot_LCRP_Q3',
    str_detect(GOterms, 'protein') & domain == 1 & LCR_Len_bin == 'Q4' ~ 'Prot_LCRP_Q4',
    str_detect(GOterms, 'protein') & domain == 0 & LCR_Len_bin == '0' ~ 'Prot_LCRA_0',
    str_detect(GOterms, 'protein') & domain == 0 & LCR_Len_bin == 'Q1' ~ 'Prot_LCRA_Q1',
    str_detect(GOterms, 'protein') & domain == 0 & LCR_Len_bin == 'Q2' ~ 'Prot_LCRA_Q2',
    str_detect(GOterms, 'protein') & domain == 0 & LCR_Len_bin == 'Q3' ~ 'Prot_LCRA_Q3',
    str_detect(GOterms, 'protein') & domain == 0 & LCR_Len_bin == 'Q4' ~ 'Prot_LCRA_Q4'),
    
    Lipid = case_when(
    str_detect(GOterms, 'lipid') & domain == 1 & LCR_Len_bin == '0' ~ 'Lip_LCRP_0',
    str_detect(GOterms, 'lipid') & domain == 1 & LCR_Len_bin == 'Q1' ~ 'Lip_LCRP_Q1',
    str_detect(GOterms, 'lipid') & domain == 1 & LCR_Len_bin == 'Q2' ~ 'Lip_LCRP_Q2',
    str_detect(GOterms, 'lipid') & domain == 1 & LCR_Len_bin == 'Q3' ~ 'Lip_LCRP_Q3',
    str_detect(GOterms, 'lipid') & domain == 1 & LCR_Len_bin == 'Q4' ~ 'Lip_LCRP_Q4',
    str_detect(GOterms, 'lipid') & domain == 0 & LCR_Len_bin == '0' ~ 'Lip_LCRA_0',
    str_detect(GOterms, 'lipid') & domain == 0 & LCR_Len_bin == 'Q1' ~ 'Lip_LCRA_Q1',
    str_detect(GOterms, 'lipid') & domain == 0 & LCR_Len_bin == 'Q2' ~ 'Lip_LCRA_Q2',
    str_detect(GOterms, 'lipid') & domain == 0 & TM_Len_bin == 'Q3' ~ 'Lip_LCRA_Q3',
    str_detect(GOterms, 'lipid') & domain == 0 & TM_Len_bin == 'Q4' ~ 'Lip_LCRA_Q4')) %>% 
  select('prot_id', 'tAI', 'RNA', 'DNA', 'Mito', 'Prot', 'Lipid', 'LCR_length') %>%  # removed extra comma here
  pivot_longer(cols=c('RNA', 'DNA', 'Mito', 'Prot', 'Lipid'),
               names_to = 'GO_term',
               values_to = 'Presence') %>% 
  filter(!is.na(Presence))




PFAM_LCR_GOterms <- chosen_GOterms %>% 
  mutate(
    RNA = case_when(
      str_detect(GOterms, 'RNA') & str_detect(PFAM_LCR, '.') ~ 'RNA_PFAM_LCRP',
      TRUE ~ 'RNA_PFAM_LCRA'),
    DNA = case_when(
      str_detect(GOterms, 'DNA') & str_detect(PFAM_LCR, '.') ~ 'DNA_PFAM_LCRP',
      TRUE ~ 'DNA_PFAM_LCRA'),
    Mito = case_when(
      str_detect(GOterms, 'mitochondrial') & str_detect(PFAM_LCR, '.') ~ 'Mito_PFAM_LCRP',
      TRUE ~ 'Mito_PFAM_LCRA'),
    Prot = case_when(
      str_detect(GOterms, 'protein') & str_detect(PFAM_LCR, '.') ~ 'Prot_PFAM_LCRP',
      TRUE ~ 'Prot_PFAM_LCRA'),
    Lipid = case_when(
      str_detect(GOterms, 'lipid') & str_detect(PFAM_LCR, '.') ~ 'Lip_PFAM_LCRP',
      TRUE ~ 'Lip_PFAM_LCRA'),
    ) %>% 
  select('prot_id', 'tAI', 'RNA', 'DNA', 'Mito', 'Prot', 'Lipid') %>% 
  pivot_longer(cols=c('RNA', 'DNA', 'Mito', 'Prot', 'Lipid'),
               names_to = 'GO_term',
               values_to = 'Presence')


# phylum


#```


#Zależność między obecnością domen w białku a wartością tAI tego białka. Grupuje białka ze względu na obecność domen.
#```{r domain presence violin}

colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(domain_GOterms$Presence)))

equation_label <- create_statistical_label(domain_GOterms, domain_GOterms$tAI, domain_GOterms$Presence, TRUE)

count_values <- domain_GOterms %>% 
    group_by(Presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(Presence, "\nn =", count)) %>% 
    pull(label)



GOterms_domain_boxplot <- ggplot(domain_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
  geom_boxplot(width = 0.5, alpha = 0.8, notch = TRUE, color = 'black', linewidth = 0.5)+
  labs(subtitle = equation_label, x = '')+
  theme_stata()+
  theme(legend.position = 'none',
        plot.subtitle = element_text(hjust=0),
        axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
        axis.text.y=element_text(face='bold'),
        panel.grid.major.y = element_line(color='gray', size=0.5),
        panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
  scale_fill_manual(values=colors)+
  scale_x_discrete(labels=x_labels)+
  scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))



GOterms_domain_violin <- ggplot(domain_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black', size=0.5)+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = colors)+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))
  
ggsave(paste0(output_path, 'GOterms_domain_boxplot.png'), plot=GOterms_domain_boxplot, width = 12, height = 6, dpi=300)
ggsave(paste0(output_path, 'GOterms_domain_violin.png'), plot=GOterms_domain_violin, width = 12, height = 6, dpi=300)

#Zgadza sie, jest to prawda. To nie jest najlepszy wykres, wyjaśnię to opisując powstawanie tych grup.
#DNA domain presence powstało jeśli biąłko ma domenę oraz w polu 'go terms' znajduje się string 'DNA'
#Jeśli białko ma domenę, a w polu 'go terms' nie ma stringu 'DNA' to jest mu przypisana etykieta DNA_doma (absence)
#Wszystkie białka w odrębie danego gotermsu sumują się do 319200 czyli do łącznej sumy.
#Poprawne etykiety powinny być:
#DNA_goterms_presence, DNA_goterms_absence
#Czyli jeszcze raz, DNA_goterms_presence zbiera te białka, które posiadają goterm ze słowem 'DNA'. W DNA_goterms_absence sa te gotermsy, które nie spełniają tego warunku (ale mogą mieć inne gotermsy)
#Pozostaje mi zrobić inny podział w kolumnach. Muszę uwzględnić białka, które mają min 2. różńe gotermsy i niech wyświetlają się na wykresie jako osobne kolumny. Poprawi to czytelność i wiarygodność danych 
#```

#Zależnosć liczby regionów LCR w białku w zależności od wartości tAI tego białka. Tym białkom nadałem etykiety: LCR absence, Jeden LCR oraz więcej niż jeden LCR
##Stworzyć dane, w których badam samo tAI tych odcinków LCR. Muszę przekodować sekwencje białkowe LCR na nukleotydowe a następnie obliczyć z nich tAI i porównać do poziomu tAI białka
#```{r number of LCR}

colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(LCR_GOterms$Presence)))

equation_label <- create_statistical_label(LCR_GOterms, LCR_GOterms$tAI, LCR_GOterms$Presence, TRUE)

count_values <- LCR_GOterms %>% 
    group_by(Presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(Presence, "\nn =", count)) %>% 
    pull(label)



GOterms_LCRpres_boxplot <- ggplot(LCR_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
  geom_boxplot(width = 0.5, alpha = 0.8, notch = TRUE, color = 'black', linewidth = 0.5)+
  labs(subtitle = equation_label, x = '')+
  theme_stata()+
  theme(legend.position = 'none',
        plot.subtitle = element_text(hjust=0),
        axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
        axis.text.y=element_text(face='bold'),
        panel.grid.major.y = element_line(color='gray', size=0.5),
        panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
  scale_fill_manual(values=colors)+
  scale_x_discrete(labels=x_labels)+
  scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))



GOterms_LCRpres_violin <- ggplot(LCR_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black', size=0.5)+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = colors)+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))
  
ggsave(paste0(output_path, 'GOterms_LCRpres_boxplot.png'), plot=GOterms_LCRpres_boxplot, width = 12, height = 6, dpi=300)
ggsave(paste0(output_path, 'GOterms_LCRpres_violin.png'), plot=GOterms_LCRpres_violin, width = 12, height = 6, dpi=300)

#```

#Zależnosć długości odcinka LCR w białku w wartości tAI tego białka. Tym białkom nadałem etykiety: dłguości odcinków LCR pogrupowałem w 4 biny zgodnie z kwartylami
##Stworzyć dane, w których badam samo tAI tych odcinków LCR. Muszę przekodować sekwencje białkowe LCR na nukleotydowe a następnie obliczyć z nich tAI i porównać do poziomu tAI białka
#```{r vs LCR presence}
colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(LCR_len_GOterms$Presence)))

equation_label <- create_statistical_label(LCR_len_GOterms, LCR_len_GOterms$tAI, LCR_len_GOterms$Presence, TRUE)

count_values <- LCR_len_GOterms %>% 
    group_by(Presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(Presence, "\nn =", count)) %>% 
    pull(label)



GOterms_LCRlen_boxplot <- ggplot(LCR_len_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
  geom_boxplot(width = 0.5, alpha = 0.8, notch = TRUE, color = 'black', linewidth = 0.5)+
  labs(subtitle = equation_label, x = '')+
  theme_stata()+
  theme(legend.position = 'none',
        plot.subtitle = element_text(hjust=0),
        axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
        axis.text.y=element_text(face='bold'),
        panel.grid.major.y = element_line(color='gray', size=0.5),
        panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
  scale_fill_manual(values=colors)+
  scale_x_discrete(labels=x_labels)+
  scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))



GOterms_LCRlen_violin <- ggplot(LCR_len_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black', size=0.5)+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = colors)+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))

ggsave(paste0(output_path, 'GOterms_LCRlen_boxplot.png'), plot=GOterms_LCRlen_boxplot, width = 12, height = 6, dpi=300)
ggsave(paste0(output_path, 'GOterms_LCRlen_violin.png'), plot=GOterms_LCRlen_violin, width = 12, height = 6, dpi=300)
#```

#dlugosc odcinka LCR w zaleznosci od tAI biąłka z tym LCR
#```{r LCR length scatter}
colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(LCR_len_GOterms$LCR_length)))

equation_label <- create_statistical_label(LCR_len_GOterms, LCR_len_GOterms$tAI, LCR_len_GOterms$LCR_length, TRUE)

count_values <- LCR_len_GOterms %>% 
    group_by(LCR_length) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(LCR_length, "\nn =", count)) %>% 
    pull(label)



GOterms_LCRlen_scatter <- ggplot(LCR_len_GOterms, aes(LCR_length, tAI, fill=as.factor(Presence)))+
  geom_point(aes(colour = as.factor(Presence)))+
  geom_smooth(aes(colour = as.factor(Presence)), method=loess, se = FALSE)+
  labs(subtitle = equation_label, x = 'LCR length.')+
  theme_stata()+
  theme(legend.position = 'bottom',
        plot.subtitle = element_text(hjust=0, vjust=-2, margin = margin(b = 20)),
        axis.text.x = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold'))+
  scale_fill_manual(values=colors)+
  scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))

ggsave(paste0(output_path, 'GOterms_LCRlen_scatter.png'), plot=GOterms_LCRlen_scatter, width = 12, height = 6, dpi=300)
#```

#```{r LCR length scatter}
LCR_len_GOterms_trimmed <- LCR_len_GOterms %>% 
  filter(LCR_length < 0.85 * max(LCR_length))

colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(LCR_len_GOterms_trimmed$LCR_length)))

equation_label <- create_statistical_label(LCR_len_GOterms_trimmed, LCR_len_GOterms_trimmed$tAI, LCR_len_GOterms_trimmed$LCR_length, TRUE)

count_values <- LCR_len_GOterms_trimmed %>% 
    group_by(LCR_length) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(LCR_length, "\nn =", count)) %>% 
    pull(label)



GOterms_LCRlen_TRIMMED_scatter <- ggplot(LCR_len_GOterms_trimmed, aes(LCR_length, tAI, fill=as.factor(Presence)))+
  geom_point(aes(colour = as.factor(Presence)), alpha = 0.1)+
  geom_smooth(aes(colour = as.factor(Presence)), method=loess, se = FALSE)+
  labs(subtitle = equation_label, x = 'LCR length 30% the longest outliers were rejected')+
  theme_stata()+
  theme(legend.position = 'bottom',
        plot.subtitle = element_text(hjust=0, vjust=-2, margin = margin(b = 20)),
        axis.text.x = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold'))+
  scale_fill_manual(values=colors)+
  scale_y_continuous(breaks=seq(0,1, by=0.1),
                     minor_breaks=seq(0,1, by=0.05))

ggsave(paste0(output_path, 'GOterms_LCRlen_TRIMMED_scatter.png'), plot=GOterms_LCRlen_TRIMMED_scatter, width = 12, height = 6, dpi=300)
#```

# Obecność peptydów sygnałowych w białku vs tAI tego białka
#```{r signal presence violin and hist }
colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(signal_GOterms$Presence)))

equation_label <- create_statistical_label(signal_GOterms, signal_GOterms$tAI, signal_GOterms$Presence, TRUE)

count_values <- signal_GOterms %>% 
    group_by(Presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(Presence, "\nn =", count)) %>% 
    pull(label)



GOterms_signal_presence_bosplot <- ggplot(signal_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
  geom_boxplot(width = 0.5, alpha = 0.8, notch = TRUE, color = 'black', linewidth = 0.5)+
  labs(subtitle = equation_label, x = '')+
  theme_stata()+
  theme(legend.position = 'none',
        plot.subtitle = element_text(hjust=0),
        axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
        axis.text.y=element_text(face='bold'),
        panel.grid.major.y = element_line(color='gray', size=0.5),
        panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
  scale_fill_manual(values=colors)+
  scale_x_discrete(labels=x_labels)+
  scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))



GOterms_signal_presence_violin <- ggplot(signal_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black', size=0.5)+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = colors)+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))

ggsave(paste0(output_path, 'GOterms_signal_presence_bosplot.png'), plot=GOterms_signal_presence_bosplot, width = 12, height = 6, dpi=300)
ggsave(paste0(output_path, 'GOterms_signal_presence_violin.png'), plot=GOterms_signal_presence_violin, width = 12, height = 6, dpi=300)

#```


# tAI białek vs obecność odcinków LCR zachodzących (lub nie) na domenę PFAM
#```{r PFAM domain overlaping with LCR or not overlaping violin hist}
colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(PFAM_LCR_GOterms$Presence)))

equation_label <- create_statistical_label(PFAM_LCR_GOterms, PFAM_LCR_GOterms$tAI, PFAM_LCR_GOterms$Presence, TRUE)

count_values <- PFAM_LCR_GOterms %>% 
    group_by(Presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(Presence, "\nn =", count)) %>% 
    pull(label)



GOterms_PFAM_LCR_presence_boxplot <- ggplot(PFAM_LCR_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
  geom_boxplot(width = 0.5, alpha = 0.8, notch = TRUE, color = 'black', linewidth = 0.5)+
  labs(subtitle = equation_label, x = '')+
  theme_stata()+
  theme(legend.position = 'none',
        plot.subtitle = element_text(hjust=0),
        axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
        axis.text.y=element_text(face='bold'),
        panel.grid.major.y = element_line(color='gray', size=0.5),
        panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
  scale_fill_manual(values=colors)+
  scale_x_discrete(labels=x_labels)+
  scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))



GOterms_PFAM_LCR_presence_violin <- ggplot(PFAM_LCR_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black', size=0.5)+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = colors)+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))
ggsave(paste0(output_path, 'GOterms_PFAM_LCR_presence_boxplot.png'), plot=GOterms_PFAM_LCR_presence_boxplot, width = 12, height = 6, dpi=300)
ggsave(paste0(output_path, 'GOterms_PFAM_LCR_presence_violin.png'), plot=GOterms_PFAM_LCR_presence_violin, width = 12, height = 6, dpi=300)
#```

#tAI białek z podziałem na obecność lub nieobecność elementów transmembranowych w białku
#```{r transmembrane elements presence in protein violin hist}
colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(TM_GOterms$Presence)))

equation_label <- create_statistical_label(TM_GOterms, TM_GOterms$tAI, TM_GOterms$Presence, TRUE)

count_values <- TM_GOterms %>% 
    group_by(Presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(Presence, "\nn =", count)) %>% 
    pull(label)



GOterms_TM_presence_boxplot <- ggplot(TM_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
  geom_boxplot(width = 0.5, alpha = 0.8, notch = TRUE, color = 'black', linewidth = 0.5)+
  labs(subtitle = equation_label, x = '')+
  theme_stata()+
  theme(legend.position = 'none',
        plot.subtitle = element_text(hjust=0),
        axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
        axis.text.y=element_text(face='bold'),
        panel.grid.major.y = element_line(color='gray', size=0.5),
        panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
  scale_fill_manual(values=colors)+
  scale_x_discrete(labels=x_labels)+
  scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))



GOterms_TM_presence_violin <- ggplot(TM_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black', size=0.5)+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = colors)+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))

ggsave(paste0(output_path, 'GOterms_TM_presence_boxplot.png'), plot=GOterms_TM_presence_boxplot, width = 12, height = 6, dpi=300)
ggsave(paste0(output_path, 'GOterms_TM_presence_violin.png'), plot=GOterms_TM_presence_violin, width = 12, height = 6, dpi=300)
#```


#Długość odcinka transmembranowego (z podziałem na kwartyle) vs tAI białka 
#```{r TM len bin boxplot}
colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(TM_len_GOterms$Presence)))

equation_label <- create_statistical_label(TM_len_GOterms, TM_len_GOterms$tAI, TM_len_GOterms$Presence, TRUE)

count_values <- TM_len_GOterms %>% 
    group_by(Presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(Presence, "\nn =", count)) %>% 
    pull(label)



GOterms_TM_len_boxplot <- ggplot(TM_len_GOterms, aes(Presence, tAI, fill=as.factor(Presence)))+
  geom_boxplot(width = 0.5, alpha = 0.8, notch = TRUE, color = 'black', linewidth = 0.5)+
  labs(subtitle = equation_label, x = '')+
  theme_stata()+
  theme(legend.position = 'none',
        plot.subtitle = element_text(hjust=0),
        axis.text.x = element_text(face='bold', angle = 45, hjust = 1),
        axis.text.y=element_text(face='bold'),
        panel.grid.major.y = element_line(color='gray', size=0.5),
        panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
  scale_fill_manual(values=colors)+
  scale_x_discrete(labels=x_labels)+
  scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))

ggsave(paste0(output_path, 'GOterms_TM_len_boxplot.png'), plot=GOterms_TM_len_boxplot, width = 12, height = 6, dpi=300)
#```





