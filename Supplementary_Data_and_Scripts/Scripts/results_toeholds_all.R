####################################
####    Results toeholds all    ####
#### This script makes figures  ####
#### 4 and 5 for the manuscript ####
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
library(agricolae)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

theme_set(theme_cowplot() + 
            theme(panel.background = element_rect(fill = 'white'),
                  plot.background = element_rect(fill = 'white'), 
                  axis.text = element_text(size = 18), 
                  axis.title = element_text(size = 20, face = 'bold'), 
                  strip.text = element_text(size = 20, face = 'bold'), 
                  strip.background = element_blank())
) 

## Set path to the main folder
# main_folder <- '/path/to/Supplementary_Data_and_Scripts/Data'
main_folder <- '/media/angelfcc/Angel_backup/iGEM/Toeholder_SuppData_Github/Toeholder_Data_and_Scripts/Supplementary_Data_and_Scripts/Data/'
setwd(main_folder)

#### Figure 4 ####

# Load the data
data_green_mfe <- read_delim('toehold_Green_results.tsv', delim = '\t')
data_green_on_off <- read_delim('all_toeholds_Green_2014_formatted.csv', delim = ',')

# Merge the tables
merged_data_green <- inner_join(x = data_green_mfe, 
                                y = data_green_on_off %>% 
                                  unite(Table, `Toehold switch number`, col = 'ID', sep = '_'), 
                                by = c('ID'))

merged_data_green %<>% 
  mutate(diff_bound_unbound = Toehold_trigger_MFE - Toehold_MFE - Trigger_MFE)

# Load the data
green_raw_data <- read_delim('all_toeholds_Green_ON_OFF_ONOFF_data.csv', delim = '\t', 
                             locale = locale(decimal_mark = ','))

# Load the GC data
gc_data <- read_delim('toehold_Green_gc_results.tsv', delim = '\t')

# Merge with the ON / OFF data
green_gc_merged_data <- inner_join(x = gc_data %>%
                                     separate(col = ID, into = c('Table', 'toehold'), sep = '_') %>%
                                     mutate(Table = ifelse(Table == 'S1A', 'S1B', Table)) %>%
                                     unite(Table, toehold, col = 'ID', sep = '_'), 
                                   y = green_raw_data %>% 
                                     # filter(Table != 'S3A') %>%
                                     unite(Table, `Toehold switch number`, col = 'ID', sep = '_'), 
                                   by = c('ID')) %>%
  mutate(diff_bound_unbound = Toehold_trigger_MFE - Toehold_MFE - Trigger_MFE)

data_on_off_gc_key <- green_gc_merged_data %>% 
  mutate(
    GC_strong = factor(GC_strong, levels = c(0, 1, 2, 3, 4, 5)),
    GC_weak = factor(GC_weak, levels = c(0, 1, 2, 3, 4)),
    FE_check = (substr(x = ID, start = 1, stop = 3) != 'S1B')) %>%
  pivot_longer(cols = c('ON median avg', 'OFF median avg', 'ON/OFF median avg'), 
               names_to = 'Metric', values_to = 'Metric_value') %>%
  filter(FE_check != 1) %>%
  filter(Metric == 'ON/OFF median avg')

## Do ANOVA and post hoc comparisons
m <- aov(log10(Metric_value) ~ GC_strong, data=data_on_off_gc_key)
anova_test <- anova(m)
tukey_test <- HSD.test(m, trt = 'GC_strong')
tukey_test2 <- TukeyHSD(x = m, 'GC_strong', conf.level = 0.95)

tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = c(3, 3, 3, 3, 3))
tukey_annotation$GC_strong <- rownames(tukey_annotation)

toehold_counts <- data_on_off_gc_key %>%
  ungroup() %>% group_by(GC_strong) %>%
  summarise(toehold_count = str_c('n = ', n()))

p_strong <- data_on_off_gc_key %>%
  ggplot(aes(x = GC_strong, y = log10(Metric_value), fill = Metric)) +
  geom_jitter(aes(colour = Metric, shape = FE_check), width = 0.1) +
  geom_boxplot(width = 0.75, alpha = 0.3) +
  geom_violin(width = 0.75, alpha = 0.3, scale = 'width') +
  theme(legend.position = 'none', legend.justification = 'center', 
  ) +
  geom_text(data = tukey_annotation,
            inherit.aes = F,
            aes(x = GC_strong, y = y, label = groups), 
            size = 5) +
  geom_text(data = toehold_counts, 
            inherit.aes = F, 
            aes(x = GC_strong, y = 3.2, label = toehold_count), 
            size = 5) +
  xlab('GC at most stable positions') + ylab('log10(ON / OFF median avg)') +
  scale_shape(guide = 'none') +
  scale_fill_manual(values = c('#619CFF')) +
  scale_colour_manual(values = c('#619CFF'))
p_strong

### Repeat for the GC_weak
## Do ANOVA and post hoc comparisons
m <- aov(log10(Metric_value) ~ GC_weak, data=data_on_off_gc_key)
anova_test <- anova(m)
tukey_test <- HSD.test(m, trt = 'GC_weak')
tukey_test2 <- TukeyHSD(x = m, 'GC_weak', conf.level = 0.95)


tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = c(3, 3, 3, 3, 3))
tukey_annotation$GC_weak <- rownames(tukey_annotation)

toehold_counts <- data_on_off_gc_key %>%
  ungroup() %>% group_by(GC_weak) %>%
  summarise(toehold_count = str_c('n = ', n()))

p_weak <- data_on_off_gc_key %>% 
  mutate(GC_strong = factor(GC_strong, levels = c(0, 1, 2, 3, 4, 5)),
         GC_weak = factor(GC_weak, levels = c(0, 1, 2, 3, 4)), 
         FE_check = (substr(x = ID, start = 1, stop = 3) != 'S1B')) %>%
  # pivot_longer(# cols = c('ON median avg', 'OFF median avg', 'ON/OFF median avg'), 
  #   cols = c('ON/OFF median avg'), 
  #              names_to = 'Metric', values_to = 'Metric_value') %>%
  filter(Metric == 'ON/OFF median avg') %>%
  ggplot(aes(x = GC_weak, y = log10(Metric_value), fill = Metric)) +
  geom_jitter(aes(colour = Metric, shape = FE_check), width = 0.1) +
  geom_boxplot(width = 0.75, alpha = 0.3) +
  geom_violin(width = 0.75, alpha = 0.3, scale = 'width') +
  theme(legend.position = 'none', legend.justification = 'center', 
  ) +
  geom_text(data = tukey_annotation,
            inherit.aes = F,
            aes(x = GC_weak, y = y, label = groups), 
            size = 5) +
  geom_text(data = toehold_counts, 
            inherit.aes = F, 
            aes(x = GC_weak, y = 3.2, label = toehold_count), 
            size = 5) +
  xlab('GC at least stable positions') + ylab('log10(ON / OFF median avg)') +
  scale_shape(guide = 'none') +
  scale_fill_manual(values = c('#619CFF')) +
  scale_colour_manual(values = c('#619CFF'))
p_weak


data_summary <- # data_on_off_gc_key %>% 
  green_gc_merged_data %>%
  ungroup() %>% 
  separate(col = 'ID', into = c('table', 'tmp')) %>%
  filter(table != 'S3A') %>%
  group_by(GC_weak, GC_strong) %>%
  # summarise(mean_ON_OFF = mean(Metric_value), num_toeholds = n())
  summarise(mean_ON_OFF = mean(`ON/OFF median avg`), num_toeholds = n())

#### Draw a heatmap with ComplexHeatmap ####

## Convert to a wide matrix
data_fig_4c <- data_summary %>% select(-num_toeholds) %>% ungroup() %>%
  pivot_wider(names_from = GC_strong, values_from = mean_ON_OFF)

## Move the first column to rownames
needed_rownames <- data_fig_4c$GC_weak
data_fig_4c %<>% select(-GC_weak)

## Reorder rows of matrices
data_fig_4c <- data_fig_4c[c(5, 4, 3, 2, 1),]
needed_rownames <- needed_rownames[c(5, 4, 3, 2, 1)]
rownames(data_fig_4c) <- needed_rownames

## Organize the numbers of toeholds in a similar way
num_toeholds_4c <- data_summary %>% select(-mean_ON_OFF) %>% ungroup() %>%
  pivot_wider(names_from = GC_strong, values_from = num_toeholds)

needed_rownames <- num_toeholds_4c$GC_weak
num_toeholds_4c %<>% select(-GC_weak)

## Reorder
num_toeholds_4c <- num_toeholds_4c[c(5, 4, 3, 2, 1),]
needed_rownames <- needed_rownames[c(5, 4, 3, 2, 1)]
rownames(num_toeholds_4c) <- needed_rownames

## Draw the heatmap
## Draw the heatmap
row_title_size = 20
column_title_size = 20
row_name_size = 18
column_name_size = 18
legend_title_size = 18
legend_label_size = 16
heatmap_width = 10
heatmap_height = 10

p <- Heatmap(
  as.matrix(data_fig_4c), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(0, 40, 80, 120, 160),
    colors = c('#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')),
  show_column_names = T, column_names_side = 'bottom',
  show_row_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = "GC at least stable positions",
  column_title = 'GC at most stable positions',
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=column_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 0, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=column_name_size,fontface='bold'),
  column_names_rot = 0,
  ## Add the number of toeholds
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if(is.na(num_toeholds_4c[i, j])){
      grid.text(
        'n = 0',
        x, y,
        gp = gpar(fontsize = 16, fontface = 'bold'))
    }else{
      grid.text(
        str_c('n = ', num_toeholds_4c[i, j], sep = ''),
        x, y,
        gp = gpar(fontsize = 16, fontface = 'bold'))
    }
  },
  show_heatmap_legend = T,
  heatmap_legend_param = list(
    at = c(0, 40, 80, 120, 160),
    title = "Mean ON/OFF ratio", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  )
)
p
# Put the figures together
p_fig4 <- plot_grid(p_strong, p_weak, grid.grabExpr(draw(p)),
                    labels = c('A', 'B', 'C'), label_size = 30, 
                    label_fontface = 'bold', nrow = 3)
ggsave(plot = p_fig4, device = cairo_pdf, width = 7, height = 17, dpi = 300, 
       filename = file.path(main_folder, '..', '..', 'Figures', 'Fig4_new.pdf'))
ggsave(plot = p_fig4, device = 'png', width = 7, height = 17, dpi = 300, 
       filename = file.path(main_folder, '..', '..', 'Figures', 'Fig4_new.png'))

#### Figure 5 ####

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
  
  # Get the number of total candidate toeholds
  subfolders <- list.dirs(path = input_folder, full.names = F, recursive = F)
  subfolders <- subfolders[!(subfolders == 'Figures')]
  total_toeholds <- max(as.numeric(subfolders)) + 1 # Adjust for zero-based
  
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

fig5A <- all_data_toeholds %>% 
  ggplot(aes(x = target_name, y = diff_energy)) + 
  geom_violin(fill = 'grey') +
  geom_boxplot(outlier.shape = NA, fill = 'grey', alpha = 0.5) + 
  geom_jitter(width = 0.15, alpha = 0.5) +
  ylab('Free energy difference (kcal/mol)') + xlab('Target name') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
fig5A

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
fig5B <- final_data_all %>% rowwise() %>%
  mutate(Category = factor(Category, levels = c('Stop codons', 'Partial match', 'Perfect match')),
         Percentage_label = paste(toString(round(Percentage, 1)), '%', sep = '')) %>%
  ggplot(aes(x = target_name, y = Percentage, fill = Category)) +
  scale_fill_manual(values = c('#ff4d4d', '#ffad33', '#009900')) +
  geom_bar(stat = 'identity', position = 'stack') +
  geom_text(aes(label = Percentage_label), size = 6, position = position_stack(vjust = 0.5),
            fontface = 'bold') +
  labs(fill = '') +
  theme(legend.position = 'top', legend.justification = 'center', 
        axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
  xlab('Target name') + ylab('Percentage (%)')
fig5B

# Put the two figures together
figure5 <- plot_grid(fig5A, fig5B, nrow = 2, labels = c('A', 'B'), label_size = 30)

ggsave(plot = figure5, device = cairo_pdf, width = 10, height = 14, dpi = 300, 
       filename = file.path(main_folder, '..', '..', 'Figures', 'Fig5.pdf'))
