#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This script saves some input variables needed by the toehold_master.py script

# The path to the sequence for which toeholds will be designed
input_seq = '/media/axelle/Angel_backup/Dropbox/Hiver2019/iGEM/Viral_toeholds/PR772_ORF.fasta'

# The length of the unpaired region of the recognition sequence (upstream of the hairpin)
length_unpaired = 14

# The length of the paired region of the recognition sequence (in the hairpin)
length_paired = 16

# The path to the output folder
output_folder = '/media/axelle/Angel_backup/Dropbox/Hiver2019/iGEM/Viral_toeholds/PR772_ORF_toeholds_2019-09-16'

# The reporter sequence to be simulated at the end of each toehold switch
reporter = 'AUGCGUAAA'

# The molecule type that is provided as input (RNA or DNA)
mol_type = 'DNA'

# The list of reference genomes to which candidate toeholds should be mapped
reference_list = '/media/axelle/Angel_backup/iGEM/Genomes/Genome_list.txt'

# The lower threshold of sequence identity for the hits to be retained
pct_ident = 90

# The upper threshold for the evalue of hits to be retained
evalue = 1e-6

#### RNA input specific paramters ####

# The minimum number of unpaired residues in the secondary structure of the target RNA
min_unpaired = 10
