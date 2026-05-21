## ----TODO---------------------------------------------------------------------



## ----setup, include=FALSE-----------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE)


## ----libraries----------------------------------------------------------------
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


## ----read data - local paths--------------------------------------------------
# main_path <- '/home/norbert/IBB/skrypty/tAI/Main_Dataset_splitted/'
# domain_bolean <- read.csv(paste0(main_path, 'domain_bolean_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'domain', 'tAI'))
# 
# signal_bolean <- read.csv(paste0(main_path, 'signal_bolean_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'signal', 'tAI'))
# 
# protein_length <- read.csv(paste0(main_path, 'protein_length_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'length', 'tAI'))
# 
# LCR_number <- read.csv(paste0(main_path, 'LCR_number_columns.tsv'), sep='\t', header=FALSE, col.names = c('prot_id','LCR_count', 'tAI'))
# 
# TM_count <- read.csv(paste0(main_path, 'transmembrane_count_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'TM_count_col', 'tAI'))
# 
# TM_length <- read.csv(paste0(main_path, 'total_transmembrane_length_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'TM_length_col', 'tAI'))
# 
# GOterms <- read.csv(paste0(main_path, 'GOterms_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'GOterms', 'tAI'))
# 
# PFAM_LCR <- read.csv(paste0(main_path, 'PFAM_domains_LCR_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'PFAM_LCR_col', 'tAI'))
# 
# PFAM_proteins <- read.csv(paste0(main_path, 'PFAM_domains_proteins_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'PFAM_prot', 'tAI'))
# 
# organism_acc <- read.csv(paste0(main_path, 'organism_acc_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'assembly_id', 'tAI'))

# output_path<-'/home/norbert/IBB/skrypty/tAI/Main_Dataset_splitted/'


## ----read data - mycelia paths------------------------------------------------
main_path <- '/home/norbert_s/tAI/all_chunks/Main_Dataset_splitted/data/'
domain_bolean <- read.csv(paste0(main_path, 'domain_bolean_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'domain', 'tAI'))

signal_bolean <- read.csv(paste0(main_path, 'signal_bolean_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'signal', 'tAI'))

protein_length <- read.csv(paste0(main_path, 'protein_length_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'length', 'tAI'))

LCR_number <- read.csv(paste0(main_path, 'LCR_number_columns.tsv'), sep='\t', header=FALSE, col.names = c('prot_id','LCR_count', 'tAI'))

TM_count <- read.csv(paste0(main_path, 'transmembrane_count_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'TM_count_col', 'tAI'))

TM_length <- read.csv(paste0(main_path, 'total_transmembrane_length_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'TM_length_col', 'tAI'))

GOterms <- read.csv(paste0(main_path, 'GOterms_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'GOterms', 'tAI'))

PFAM_LCR <- read.csv(paste0(main_path, 'PFAM_domains_LCR_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'PFAM_LCR_col', 'tAI'))

PFAM_proteins <- read.csv(paste0(main_path, 'PFAM_domains_proteins_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'PFAM_prot', 'tAI'))

organism_acc <- read.csv(paste0(main_path, 'organism_acc_column.tsv'), sep='\t', header=FALSE, col.names = c('prot_id', 'assembly_id', 'tAI'))

output_path<-'/home/norbert_s/tAI/all_chunks/Main_Dataset_splitted/plots/'


## ----modify data--------------------------------------------------------------
#replace numerical values (booleans 0 or 1) by strings (Presence, No domain)
domain_bolean <- domain_bolean %>%
    mutate(presence = factor(case_when(
        domain == 1 ~ 'Domain presence',
        domain == 0 ~ 'Domain absence'
), levels=c('Domain absence', 'Domain presence')))


#replace numerical values (booleans 0 or 1) by strings (Presence, absence)
signal_bolean <- signal_bolean %>%
    mutate(presence = factor(case_when(
        signal == 1 ~ 'Signal presence',
        signal == 0 ~ 'Signal absence'
        ), levels=c('Signal absence', 'Signal presence')))


protein_length$length <- 
    as.numeric(protein_length$length)
breaks = seq(1, max(protein_length$length, na.rm = TRUE), length.out = 11)
protein_length$bin <- cut(
  protein_length$length,
  breaks = breaks,
  labels = 1:10,
  right = FALSE)
protein_length$bin[is.na(protein_length$bin)] <- 9

LCR_number <- LCR_number %>%
    mutate(presence = factor(case_when(
        LCR_count ==0 ~ 'LCR absence',
        LCR_count != 0 ~'LCR presence'),
        levels=c('LCR absence', 'LCR presence')))


GOterms <- GOterms[GOterms$GOterms != "",]
GOterms <- GOterms %>%
separate_rows(GOterms, sep = "\\|")

PFAM_LCR <- PFAM_LCR %>% 
  mutate(presence = factor(case_when(
    PFAM_LCR_col == '' ~ 'PFAM LCR absence',
    PFAM_LCR_col != '' ~ 'PFAM LCR presence'),
    levels=c('PFAM LCR absence', 'PFAM LCR presence')))


PFAM_proteins <- PFAM_proteins[PFAM_proteins$PFAM_prot != "",]
PFAM_proteins <- PFAM_proteins %>%
separate_rows(PFAM_prot, sep = " ")

TM_length$TM_length_col <- 
    as.numeric(TM_length$TM_length_col)

breaks = c(0, seq(1, max(TM_length$TM_length_col, na.rm = TRUE), length.out = 11))
TM_length$bin <- cut(
  TM_length$TM_length_col,
  breaks = breaks,
  labels = 0:10,
  right = FALSE
)
TM_length$bin[is.na(TM_length$bin)] <- 9


TM_count <- TM_count %>% 
  mutate(presence = factor(case_when(
    TM_count_col == 0 ~ 'Transmembrane elements absence',
    TM_count_col != 0 ~ 'Transmembrane elements presence'),
    levels = c('Transmembrane elements absence', 'Transmembrane elements presence')))


mean_tAI <- mean(domain_bolean$tAI, na.rm=TRUE)


## ----functions----------------------------------------------------------------
create_regression_equation <- function(data, formula) {
  model <- lm(formula, data=data)
  summary_model <- summary(model)
  r_sq <- round(summary_model$r.squared, 5)
  coef <- coef(model)
  intercept <- round(coef[1], 2)
  slope <- round(coef[2], 2)
  p_value <- round(summary_model$coefficients[2, 4], 5)
  equation_label <- paste("y = ", intercept, " + ", slope, " * x,", "        R² = ", r_sq, ",        p value = ", p_value, sep = "")
  return(equation_label)
}

create_statistical_label <- function(data, var1, var2, parametric_assumption = TRUE) {
  
  if (parametric_assumption) {
    ad_test <- ad.test(var1)
    
    if (ad_test$p.value < 0.05) {
      # ANOVA test
      test = 'ANOVA p-value: '
      anova_results <- summary(aov(var1 ~ var2, data = data))
      p_value <- anova_results[[1]]$'Pr(>F)'[1]
    } else {
      # Kruskal-Wallis test
      test = 'Kruskal-Wallis test p-value: '
      kw_test <- kruskal.test(var1 ~ var2, data = data)
      p_value <- kw_test$p.value
    }
  } else {
    ad_test <- ad.test(var1)
    
    if (ad_test$p.value < 0.05) {
      # Mann-Whitney test
      test = 'Mann-Whitney test p-value: '
      mann_whitney_test <- wilcox.test(var1 ~ var2, data = data)
      p_value <- mann_whitney_test$p.value
    } else {
      # T-test
      test = 'T-student test p-value: '
      t_student_result <- t.test(var1, mu = 0.5)
      p_value <- t_student_result$p.value
    }
  }
  
  # Create the equation label
  equation_label <- paste(test, round(p_value, 5), ".")
  return(equation_label)
}


## ----domain presence violin---------------------------------------------------

equation_label <- create_statistical_label(domain_bolean, domain_bolean$tAI, domain_bolean$presence, FALSE)

count_values <- domain_bolean %>% 
    group_by(presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(presence, "\nn=", count)) %>% 
    pull(label)


#violin
domains_violin <- ggplot(domain_bolean, aes(presence, tAI, fill=as.factor(presence)))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black', size=0.5)+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold'),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = c('Domain absence' = '#696969', 'Domain presence'='#C0C0C0'))+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))
ggsave('domains_violin.png', plot=domains_violin, width = 8, height = 6, dpi=300)
#hist
domains_hist <- ggplot(domain_bolean, aes(x=tAI, color=presence, fill=presence))+
  geom_histogram(bins=30, color='black',  position='dodge', alpha=0.8)+
  scale_fill_manual(values=c('Domain absence' = '#696969', 'Domain presence'='#C0C0C0'))+
  geom_density(lwd=0.5, linetype=1, alpha=0)+
  scale_color_manual(values=c('Domain absence' = '#696969', 'Domain presence'='#C0C0C0'))+
  geom_vline(data=domain_bolean, aes(xintercept = mean_tAI), linetype = 'dashed',  linewidth = 1)+
  labs(subtitle=equation_label)+
  theme_stata()+
  theme(
    legend.position = 'bottom',
    plot.subtitle = element_text(hjust=0),
    axis.text.x = element_text(face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )
ggsave('domains_hist.png', plot=domains_hist, width = 8, height = 6, dpi=300)
  



## ----LCR presence violin and hist---------------------------------------------
equation_label <- create_statistical_label(LCR_number, LCR_number$tAI, LCR_number$presence, FALSE)

count_values <- LCR_number %>% 
    group_by(presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(presence, "\nn=", count)) %>% 
    pull(label)


#violin
LCR_presence_violin <- ggplot(LCR_number, aes(presence, tAI, fill=presence))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black')+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold'),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = c('LCR absence' = '#696969', 'LCR presence'='#C0C0C0'))+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))

LCR_presence_violin
ggsave('LCR_presence_violin.png', plot=LCR_presence_violin, width = 8, height = 6, dpi=300)

#hist
LCR_presence_hist <- ggplot(LCR_number, aes(x=tAI, color=presence, fill=presence))+
  geom_histogram(bins=30, color='black',  position='dodge', alpha=0.8)+
  scale_fill_manual(values=c('LCR absence' = '#696969', 'LCR presence'='#C0C0C0'))+
  geom_density(lwd=0.5, linetype=1, alpha=0.0)+
  scale_color_manual(values=c('LCR absence' = '#696969', 'LCR presence'='#C0C0C0'))+
  geom_vline(data=LCR_number, aes(xintercept = mean_tAI), linetype = 'dashed',  linewidth = 1)+
  labs(subtitle=equation_label)+
  theme_stata()+
  theme(
    legend.position = 'bottom',
    plot.subtitle = element_text(hjust=0),
    axis.text.x = element_text(face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )
ggsave('LCR_presence_hist.png', plot=LCR_presence_hist, width = 8, height = 6, dpi=300)


## ----number of LCR------------------------------------------------------------
colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(LCR_number$LCR_count)))
LCR_number$LCR_count <- as.factor(LCR_number$LCR_count)
equation_label <- create_statistical_label(LCR_number, LCR_number$tAI, LCR_number$LCR_count, TRUE)



count_values <- LCR_number %>% 
    group_by(LCR_count) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(LCR_count, "\nn=", count)) %>% 
    pull(label)


#violin
LCR_number_violin <- ggplot(LCR_number, aes(LCR_count, tAI, fill=as.factor(LCR_count)))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, color='black')+
    labs(subtitle=equation_label)+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold'),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = colors)+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))

ggsave('LCR_number_violin.png', plot=LCR_number_violin, width = 8, height = 6, dpi=300)


#boxplot
LCR_number_hist <- ggplot(LCR_number, aes(LCR_count, tAI, fill = LCR_count))+
  geom_boxplot(alpha=0.8, notch=TRUE, color='black')+
  labs(subtitle=equation_label)+
  theme_stata()+
  theme(legend.position = 'none',
        plot.subtitle = element_text(hjust=0),
        axis.text.x = element_text(face='bold'),
        axis.text.y=element_text(face='bold'),
        panel.grid.major.y = element_line(color='gray', size=0.5),
        panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
  scale_fill_manual(values = colors)+
  scale_x_discrete(labels=x_labels)+
  scale_y_continuous(breaks=seq(0,1, by=0.1),
                     minor_breaks=seq(0,1, by=0.05))

ggsave('LCR_number_hist.png', plot=LCR_number_hist, width = 8, height = 6, dpi=300)


## ----signal presence violin and hist------------------------------------------
equation_label <- create_statistical_label(signal_bolean, signal_bolean$tAI, signal_bolean$presenc, FALSE)

count_values <- signal_bolean %>% 
    group_by(presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(presence, "\nn=", count)) %>% 
    pull(label)


#violin
signal_presence_violin <- ggplot(signal_bolean, aes(presence, tAI, fill=presence))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black')+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold'),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = c('Signal absence' = '#696969', 'Signal presence'='#C0C0C0'))+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))


ggsave('signal_presence_violin.png', plot=signal_presence_violin, width = 8, height = 6, dpi=300)
#hist
signal_presence_hist <- ggplot(signal_bolean, aes(x=tAI, color=presence, fill=presence))+
  geom_histogram(bins=30, color='black',  position='dodge', alpha=0.8)+
  scale_fill_manual(values=c('Signal absence' = '#696969', 'Signal presence'='#C0C0C0'))+
  geom_density(lwd=0.5, linetype=1, alpha=0.0)+
  scale_color_manual(values=c('Signal absence' = '#696969', 'Signal presence'='#C0C0C0'))+
  geom_vline(data=LCR_number, aes(xintercept = mean_tAI), linetype = 'dashed',  linewidth = 1)+
  labs(subtitle=equation_label)+
  theme_stata()+
  theme(
    legend.position = 'bottom',
    plot.subtitle = element_text(hjust=0),
    axis.text.x = element_text(face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )

ggsave('signal_presence_hist.png', plot=signal_presence_violin, width = 8, height = 6, dpi=300)


## ----protein length scatter---------------------------------------------------
#equation_label <- create_regression_equation(protein_length, tAI~length)


#p <- ggplot(protein_length, aes(length, tAI))+
#  geom_point(alpha=0.3)+
#  geom_smooth(method=loess, se=FALSE, color='696969')+
#  labs(subtitle=equation_label, x='Protein length')+
#  theme_stata()+
#  theme(
#    legend.position = 'bottom',
#    plot.subtitle = element_text(hjust=0),
#    axis.text.x = element_text(face = 'bold'),
#    axis.text.y = element_text(face = 'bold')
#  ) +
#  scale_y_continuous(breaks=seq(0,1, by=0.1),
#                       minor_breaks=seq(0,1, by=0.05))


#protein_length_scatter <- ggMarginal(p, type='densigram', size=5, fill='#696969')
#ggsave('protein_length_scatter.png', plot=protein_length_scatter, width = 8, height = 6, dpi=300)


## ----protein bins length boxplot----------------------------------------------
colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(protein_length$length)))
equation_label <- create_statistical_label(protein_length, protein_length$tAI, protein_length$length, TRUE)

count_values <- protein_length %>%
    group_by(bin) %>%
    summarise(count = n())


protein_length_boxplot <- ggplot(protein_length, aes(x=factor(bin), y=tAI, fill=factor(bin)))+
    geom_boxplot(alpha=0.8, outlier.shape=16, outlier.size=1, width=0.7)+
    geom_point(color='black', size=1, shape=16, alpha=0.5)+
    labs(subtitle=equation_label, x='Proteins Grouped into Bins with 10% Intervals')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x=element_text(face='bold', angle=60, hjust=1),
          axis.text.y=element_text(face='bold'))+
    scale_color_manual(values = colors)+
    scale_x_discrete(labels = function(x) {
        count_values$count[count_values$bin %in% as.integer(x) ] %>% 
            paste(x, ", n=", .)})
ggsave('protein_length_boxplot.png', plot=protein_length_boxplot, width = 8, height = 6, dpi=300)


## ----PFAM domain overlaping with LCR or not overlaping violin hist------------

equation_label <- create_statistical_label(PFAM_LCR, PFAM_LCR$tAI, PFAM_LCR$presence, FALSE)

count_values <- PFAM_LCR %>% 
    group_by(presence) %>% 
    summarise(count=n())
x_labels <- count_values %>% 
    mutate(label = paste(presence, "\nn=", count)) %>% 
    pull(label)


#violin
PFAM_LCR_violin <-ggplot(PFAM_LCR, aes(presence, tAI, fill=presence))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black')+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold'),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = c('PFAM LCR absence' = '#696969', 'PFAM LCR presence'='#C0C0C0'))+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))


ggsave('PFAM_LCR_violin.png', plot=PFAM_LCR_violin, width = 8, height = 6, dpi=300)

#hist
PFAM_LCR_hist <- ggplot(PFAM_LCR, aes(x=tAI, color=presence, fill=presence))+
  geom_histogram(bins=30, color='black',  position='dodge', alpha=0.8)+
  scale_fill_manual(values=c('PFAM LCR absence' = '#696969', 'PFAM LCR presence'='#C0C0C0'))+
  geom_density(lwd=0.5, linetype=1, alpha=0.0)+
  scale_color_manual(values=c('PFAM LCR absence' = '#696969', 'PFAM LCR presence'='#C0C0C0'))+
  geom_vline(data=LCR_number, aes(xintercept = mean_tAI), linetype = 'dashed',  linewidth = 1)+
  labs(subtitle=equation_label)+
  theme_stata()
  theme(
    legend.position = 'bottom',
    plot.subtitle = element_text(hjust=0),
    axis.text.x = element_text(face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )

ggsave('PFAM_LCR_hist.png', plot=PFAM_LCR_hist, width = 8, height = 6, dpi=300)


## ----transmembrane elements presence in protein violin hist-------------------

equation_label <- create_statistical_label(TM_count, TM_count$tAI, TM_count$presence, FALSE)

count_values <- TM_count %>% 
    group_by(presence) %>% 
    summarise(count=n())
x_labels <- count_values %>%
  mutate(label = paste(presence, "\nn=", count_values$count)) %>% 
  pull(label)


#violin
TM_count_violin <- ggplot(TM_count, aes(presence, tAI, fill=presence))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, alpha=0.8, notch=TRUE, fill='#FFFFF0', color='black')+
    labs(subtitle=equation_label, x='')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x = element_text(face='bold'),
          axis.text.y=element_text(face='bold'),
          panel.grid.major.y = element_line(color='gray', size=0.5),
          panel.grid.minor.y = element_line(color='lightgray', size=0.3))+
    scale_fill_manual(values = c('Transmembrane elements absence' = '#696969', 'Transmembrane elements presence'='#C0C0C0'))+
    scale_x_discrete(labels=x_labels)+
    scale_y_continuous(breaks=seq(0,1, by=0.1),
                       minor_breaks=seq(0,1, by=0.05))

ggsave('TM_count_violin.png', plot=TM_count_violin, width = 8, height = 6, dpi=300)

#hist
TM_count_hist <- ggplot(TM_count, aes(x=tAI, color=presence, fill=presence))+
  geom_histogram(bins=30, color='black',  position='dodge', alpha=0.8)+
  scale_fill_manual(values=c('Transmembrane elements absence' = '#696969', 'Transmembrane elements presence'='#C0C0C0'))+
  geom_density(lwd=0.5, linetype=1, alpha=0.0)+
  scale_color_manual(values=c('Transmembrane elements absence' = '#696969', 'Transmembrane elements presence'='#C0C0C0'))+
  geom_vline(data=LCR_number, aes(xintercept = mean_tAI), linetype = 'dashed',  linewidth = 1)+
  labs(subtitle=equation_label)+
  theme_stata()+
  theme(
    plot.subtitle = element_text(hjust=0),
    legend.position = 'bottom',
    axis.text.x = element_text(face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )
ggsave('TM_count_hist.png', plot=TM_count_hist, width = 8, height = 6, dpi=300)


## ----transmembrane elements in proteins boxplot-------------------------------
colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(TM_count$TM_count_col)))
equation_label <- create_statistical_label(TM_count, TM_count$tAI, TM_count$TM_count_col, TRUE)



count_transmembrane <- TM_count %>%
  group_by(TM_count_col) %>%
  summarise(count = n())

TM_count_boxplot <- ggplot(TM_count, aes(x=factor(TM_count_col), y=tAI, fill=factor(TM_count_col)))+
    geom_boxplot(alpha=0.8, outlier.shape=16, outlier.size=1, width=0.7)+
    geom_point(color='black', size=1, shape=16, alpha=0.5)+
    labs(subtitle=equation_label, x='Counts of transmembrane elements in protein')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x=element_text(face='bold', angle=60, hjust=1),
          axis.text.y=element_text(face='bold'))+
    scale_color_manual(values = colors)+
    scale_x_discrete(labels = function(x) {
        count_transmembrane$count[count_transmembrane$TM_count_col %in% as.integer(x)] %>%
            paste(x, ", n=", .)
  })
ggsave('TM_count_boxplot.png', plot=TM_count_boxplot, width = 8, height = 6, dpi=300)


## ----length of transmembrane elements in protein boxplot----------------------
colors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(TM_length$TM_length_col)))
equation_label <- create_statistical_label(TM_length, TM_length$tAI, TM_length$TM_length_col, TRUE)

count_transmembrane_bins <- TM_length %>%
    group_by(bin) %>%
    summarise(count = n())


TM_length_boxplot <- ggplot(TM_length, aes(x=factor(bin), y=tAI, fill=factor(bin)))+
    geom_boxplot(alpha=0.8, outlier.shape=16, outlier.size=1, width=0.7)+
    geom_point(color='black', size=1, shape=16, alpha=0.5)+
    labs(subtitle=equation_label, x='Length of transmembrane elements in protein')+
    theme_stata()+
    theme(legend.position = 'none',
          plot.subtitle = element_text(hjust=0),
          axis.text.x=element_text(face='bold', angle=60, hjust=1),
          axis.text.y=element_text(face='bold'))+
    scale_color_manual(values = colors)+
    scale_x_discrete(labels = function(x) {
        count_transmembrane_bins$count[count_transmembrane_bins$bin %in% as.integer(x) ] %>% 
            paste(x, ", n=", .)})
ggsave('TM_length_boxplot.png', plot=TM_length_boxplot, width = 8, height = 6, dpi=300)


## ----violin plots with every binary value-------------------------------------
y_limist<- c(0, 1)
joined <- plot_grid(domains_violin,
          LCR_presence_violin +theme(axis.text.y = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.title.y = element_blank()),
          signal_presence_violin, 
          PFAM_LCR_violin, 
          TM_count_violin, 
          nrow=2, ncol=3, align="v",
          labels='auto')

ggsave('joined.png', plot=joined, width = 8, height = 6, dpi=300)
#Usunac elementy z y osi dla pozostalych dwoch wykresow
# zmniejszyc czcionke
# zmniejszyc odstep miedzy liniami na osi y
# zamiast mann whitney napisac samo p-value
# nic nie poprawiac zapytac ani czy taka forma jest git

