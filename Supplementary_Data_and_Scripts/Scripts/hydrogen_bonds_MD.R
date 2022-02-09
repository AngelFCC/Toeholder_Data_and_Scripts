##############################################
####        Figure hydrogen bonds         ####
#### This script allows me to draw the    ####
#### figures on the counts of hydrogen    ####
#### bonds to analyze the stability of    ####
#### toeholds during molecular dynamics   ####
#### simulations.                         ####
##############################################

# Load libraries
library(tidyverse)
library(magrittr)
library(ggplot2)
library(cowplot)
library(readr)
library(magick)
library(rsvg)
theme_set(theme_cowplot())

# Set working directory to the main folder with the data
setwd('/path/to/Supplementary_Data_and_Scripts/Data/MD_toehold1')

# Load the file with the data for the hydrogen bonds in the hairpin
all_atom <- read_delim('hbonds_40ns_hairpin.dat', delim = ' ', col_names = c('Frame', 'Hbonds')) %>%
  mutate(Simulation = 'Explicit solvent') %>% mutate(Time = Frame / 100)

# Get the mean and standard deviation
summary_all_atom <- all_atom %>% summarise(mean = mean(Hbonds), std_dev = sd(Hbonds))

# Show the data as a figure
fig3A <- all_atom %>% mutate(Time = Frame / 100) %>%
  ggplot(aes(x = Time, y = Hbonds)) +
  geom_point() + geom_line() +
  annotate("rect", xmin=-Inf, xmax=Inf,
           ymin = summary_all_atom$mean - summary_all_atom$std_dev,
           ymax = summary_all_atom$mean + summary_all_atom$std_dev,
           alpha=0.5, fill="grey") +
  geom_abline(slope = 0, intercept = summary_all_atom$mean, colour = 'black', linetype = 'dashed') +
  ylim(0,40) +
  xlab('Time (ns)') + ylab('Hydrogen bond count')
fig3A

ggsave(fig3A, filename = '../../Figures/Fig3A.pdf', device = cairo_pdf, 
       width = 10, height = 7, dpi = 500)

#### Looking at the most stable hydrogen bonds throughout the simulation ####

# Read the data on the occupancy of each hydrogen bond
hbonds_details_raw <- read_delim('hbonds-details_no_pct_2019-10-10.dat', delim = '\t', skip = 1) 

hbonds_details <- hbonds_details_raw %>%
  mutate(donor = parse_number(donor), acceptor = parse_number(acceptor)) %>%
  rowwise() %>% mutate(res1 = min(donor, acceptor), res2 = max(donor, acceptor))

# Read the intended matches
hbonds_intended <- read_delim('pairlist_toehold1_hairpin.txt', delim = '\t', col_names = c('res1', 'res2'))

# Filter for those in the hairpin 
hbonds_details %<>% mutate(intended = ifelse(res1 == hbonds_intended$res1 && res2 == hbonds_intended$res2, 'Yes', 'No'))

# A function to label hbonds as intended or not
call_intended_hbonds <- function(hbonds_details, hbonds_intended){
  intended = c()
  for(row_hbond in 1:nrow(hbonds_details)){
    new_val = 'No'
    
    res1 = hbonds_details$res1[row_hbond]
    res2 = hbonds_details$res2[row_hbond]
    
    for(row_intended in 1:nrow(hbonds_intended)){
      # Check if this pair is in the list of intended pairs
      intended_res1 = hbonds_intended$res1[row_intended]
      intended_res2 = hbonds_intended$res2[row_intended]
      
      if(res1 == intended_res1 && res2 == intended_res2){
        new_val = 'Yes'
      }
      
    }
    intended = c(intended, new_val)
  }
  
  hbonds_details$intended <- intended
  return(hbonds_details)
}

hbonds_details <- call_intended_hbonds(hbonds_details, hbonds_intended)

#### Show a barplot of the occupancy ####

# A copy of the same data but listing the residue on the right side as res1
hbonds_right_side <- hbonds_details %>%
  select(donor, acceptor, occupancy, res2, res1, intended)

# Rename the columns
colnames(hbonds_right_side) <- c('donor', 'acceptor', 'occupancy', 'res1', 'res2', 'intended')

# Concatenate the positions for the left side with those on the right side
new_hbonds <- rbind(hbonds_details, hbonds_right_side)

# Calculate the mean occupancy for each A-T or G-C pair
new_hbonds_details <- new_hbonds %>% 
  group_by(res1, intended) %>%
  summarise(mean_occupancy = mean(occupancy)) %>%
  arrange(desc(mean_occupancy))

fig3B <- new_hbonds_details %>% 
  ungroup() %>%
  mutate(intended = ifelse(intended == 'Yes', 'Canonical', 'Not canonical')) %>%
  mutate(intended = factor(intended, levels = c('Not canonical', 'Canonical'))) %>%
  ggplot(aes(x = res1, y = mean_occupancy, fill = intended)) +
  labs(fill = '') +
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('Position') + ylab('Mean occupancy (%)')
fig3B
ggsave(fig3B, filename = '../../Figures/Barplot_occupancy_new.pdf', device = cairo_pdf, 
       width = 10, height = 7, dpi = 500)

