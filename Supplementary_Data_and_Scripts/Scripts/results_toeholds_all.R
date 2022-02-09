####################################
####    Results toeholds all    ####
#### This script will help me   ####
#### make some plots to check   ####
#### the toeholds and their     ####
#### success.                   ####
####################################

# Load libraries
library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(magrittr)
library(Cairo)
library(ggpubr)
library(rstatix)

## Need to point to the main folder that contains the data
main_folder <- '/path/to/Supplementary_Data_and_Scripts/Data'
setwd(main_folder)

# Load the data and prepare an output folder (probably better to just receive the folders here)
# oxyR_new
oxyR_folder <- 'oxyR_toeholds_2019-09-23'

# Phi6_CDS
phi6_CDS_1_folder <- 'Phi6_CDS_toeholds_2019-09-16'

# Phi_CDS_2
phi6_CDS_2_folder <- 'Phi6_CDS_2_2_toeholds_new'

# PR772_ORF
pr772_folder <- 'PR772_ORF_toeholds_2019-09-16'

# Norovirus
norovirus_folder <- 'Norovirus_toeholds_github'

# Measles
measles_folder <- 'Measlesvirus_toeholds_github'

# Alphaherpesvirus
alphaherpesvirus_folder <- 'Herpesvirus_toeholds_2_github'

#### A function to load the data, based on the section below ####

prepare_data <- function(input_folder, tag){
  
  toehold_results <- read_delim(file.path(input_folder, 'all_toeholds_results.txt'), delim = '\t')
  output_folder <- file.path(input_folder, 'Figures')
  
  if(!(dir.exists(output_folder))){
    dir.create(output_folder)
  }
  
  # Get the number of total candidate toeholds
  subfolders <- list.dirs(path = input_folder, full.names = F, recursive = F)
  subfolders <- subfolders[!(subfolders == 'Figures')]
  total_toeholds <- length(subfolders)
  
  # Count the number of failed toeholds
  stop_codons <- total_toeholds - nrow(toehold_results)
  
  # Remove NAs
  toehold_results_noNA <- toehold_results %>% filter(!(is.na(Binding_energy_toehold_mRNA)),
                                                     !(is.na(Percentage_correct_matches)),
                                                     !(is.na(Binding_energy_mRNA)))
  
  # Add a column with the sum of the individual binding energies of the toehold and the mRNA
  toehold_results_noNA %<>% mutate(diff_energy = Binding_energy_toehold_mRNA - (Binding_energy_toehold + Binding_energy_mRNA),
                                   target_name = tag)
  
  results <- list(toehold_results_noNA, stop_codons)
  
  return(results)
}

## Run the function with the datasets

# oxyR
oxyR_data <- prepare_data(oxyR_folder, 'oxyR')

# Phi6_CDS
phi6_CDS_1_data <- prepare_data(phi6_CDS_1_folder, 'phi6_CDS1')

# Phi_CDS_2
phi6_CDS_2_data <- prepare_data(phi6_CDS_2_folder, 'phi6_CDS2')

# PR772_ORF
pr772_data <- prepare_data(pr772_folder, 'pr772')

# Norovirus
norovirus_data <- prepare_data(norovirus_folder, 'Norovirus')

# Measles
measles_data <- prepare_data(measles_folder, 'Measles virus')

# Alphaherpesvirus
alphaherpesvirus_data <- prepare_data(alphaherpesvirus_folder, 'Herpesvirus')

#### Put all the data together and plot the free energy distributions ####

all_data_toeholds <- rbind(oxyR_data[[1]], phi6_CDS_1_data[[1]], phi6_CDS_2_data[[1]], pr772_data[[1]], 
                           norovirus_data[[1]], measles_data[[1]], alphaherpesvirus_data[[1]])

fig4A <- all_data_toeholds %>% 
  ggplot(aes(x = target_name, y = diff_energy)) + 
  geom_violin(fill = 'grey') +
  geom_boxplot(outlier.shape = NA, fill = 'grey', alpha = 0.5) + 
  geom_jitter(width = 0.15, alpha = 0.5) +
  ylab('Free energy difference (kcal/mol)') + xlab('Target name')
fig4A
ggsave(plot = fig4A, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = file.path(main_folder, '..', 'Figures', 'Fig4A.pdf'))

#### Prepare stacked barplots of the percentage of the different kinds of matches ####

column_tags <- c('Herpesvirus', 'Measles virus', 'Norovirus', 'oxyR', 'phi6_CDS1', 'phi6_CDS2', 'pr772')
column_stop_counts <- c(alphaherpesvirus_data[[2]], measles_data[[2]], norovirus_data[[2]],
                        oxyR_data[[2]], phi6_CDS_1_data[[2]], phi6_CDS_2_data[[2]], pr772_data[[2]])

data_stop <- as.data.frame(cbind(column_tags, column_stop_counts))
colnames(data_stop) <- c('target_name', 'count')

data_stop %<>% mutate(Category = 'Stop codons') %>%
  select(Category, target_name, count)

#### An alternative version of the distribution of binding success rates ####

toehold_results_new <- all_data_toeholds %>%
  mutate(Category = ifelse(Percentage_correct_matches == 100, 'Perfect match',
                           'Partial match')) %>%
  group_by(Category, target_name) %>%
  summarise(count = n())

# Add the counts of stop codons
all_toehold_results_new <- rbind(toehold_results_new %>% ungroup(), data_stop %>% ungroup()) %>%
  mutate(count = as.numeric(count))

# Get the total number of toeholds for each
all_toehold_totals <- all_toehold_results_new %>%
  group_by(target_name) %>%
  summarise(total = sum(count))

# Use an inner join to add the total counts
final_data_all <- inner_join(x = all_toehold_results_new, y = all_toehold_totals, 
                             by = c('target_name' = 'target_name'))

final_data_all %<>% mutate(Percentage = count * 100/ total)

# Code for the stacked barplot
fig4B <- final_data_all %>% rowwise() %>%
  mutate(Category = factor(Category, levels = c('Stop codons', 'Partial match', 'Perfect match')),
         Percentage_label = paste(toString(round(Percentage, 1)), '%', sep = '')) %>%
  ggplot(aes(x = target_name, y = Percentage, fill = Category)) +
  scale_fill_manual(values = c('#ff4d4d', '#ffad33', '#009900')) +
  geom_bar(stat = 'identity', position = 'stack') +
  geom_text(aes(label = Percentage_label), size = 6, position = position_stack(vjust = 0.5),
            fontface = 'bold') +
  labs(fill = '') +
  theme(legend.position = 'top', legend.justification = 'center') +
  xlab('Target name') + ylab('Percentage (%)')
fig4B
ggsave(plot = fig4B, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = file.path(main_folder, '..', 'Figures', 'Fig4B.pdf'))

# Put the two figures together
figure4 <- plot_grid(fig4A, fig4B, nrow = 2, labels = c('A', 'B'), label_size = 30)

ggsave(plot = figure4, device = cairo_pdf, width = 10, height = 14, dpi = 300, 
       filename = file.path(main_folder, '..', 'Figures', 'Figure4.pdf'))
