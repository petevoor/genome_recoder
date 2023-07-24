#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 19:26:11 2023

@author: petevoor
"""

import pandas as pd


# Load the CUTG csv file into a DataFrame
df = pd.read_csv('/Users/petevoor/Desktop/T7_genome_recoder/codon_usage.csv')

# Define a mapping from RNA codons to DNA codons
rna_to_dna = {
    'UUU':'TTT', 'UUC':'TTC', 'UUA':'TTA', 'UUG':'TTG', 'CUU':'CTT', 'CUC':'CTC', 
    'CUA':'CTA', 'CUG':'CTG', 'AUU':'ATT', 'AUC':'ATC', 'AUA':'ATA', 'AUG':'ATG', 
    'GUU':'GTT', 'GUC':'GTC', 'GUA':'GTA', 'GUG':'GTG', 'GCU':'GCT', 'GCC':'GCC', 
    'GCA':'GCA', 'GCG':'GCG', 'CCU':'CCT', 'CCC':'CCC', 'CCA':'CCA', 'CCG':'CCG', 
    'UGG':'TGG', 'GGU':'GGT', 'GGC':'GGC', 'GGA':'GGA', 'GGG':'GGG', 'UCU':'TCT', 
    'UCC':'TCC', 'UCA':'TCA', 'UCG':'TCG', 'AGU':'AGT', 'AGC':'AGC', 'ACU':'ACT', 
    'ACC':'ACC', 'ACA':'ACA', 'ACG':'ACG', 'UAU':'TAT', 'UAC':'TAC', 'CAA':'CAA', 
    'CAG':'CAG', 'AAU':'AAT', 'AAC':'AAC', 'UGU':'TGT', 'UGC':'TGC', 'CAU':'CAT', 
    'CAC':'CAC', 'AAA':'AAA', 'AAG':'AAG', 'CGU':'CGT', 'CGC':'CGC', 'CGA':'CGA', 
    'CGG':'CGG', 'AGA':'AGA', 'AGG':'AGG', 'GAU':'GAT', 'GAC':'GAC', 'GAA':'GAA', 
    'GAG':'GAG', 'UAA':'TAA', 'UAG':'TAG', 'UGA':'TGA'
}

# Apply the mapping to the column headers
df.rename(columns=rna_to_dna, inplace=True)

print(df)

# Write DataFrame to .csv file
df.to_csv('dna_codon_usage.csv', index=False)
