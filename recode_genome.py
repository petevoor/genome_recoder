#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 22:18:05 2023

@author: petevoor
"""

import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
import Bio.Data.CodonTable
from Bio.SeqRecord import SeqRecord
from dnachisel import DnaOptimizationProblem, AvoidPattern, EnforceTranslation, Location
from Bio.Restriction import *
import re

# Load the CUTG csv file into a DataFrame
df = pd.read_csv('dbs/dna_codon_usage.csv')

def get_codon_usage(species):
    # Query the DataFrame for the given species
    df_species = df[df['SpeciesName'] == species]
    #df_species.to_csv('df_species.csv', index=False)
    codon_usage = df_species.iloc[0, 5:].to_dict()
    
    # Convert values to float
    codon_usage = {k: float(v) for k, v in codon_usage.items()}
    
    # If the species is not found in the DataFrame, return None
    if df_species.empty:
        return None

    # Convert the DataFrame row into a dictionary and return it
    return codon_usage

def compute_relative_codon_usage(codon_usage, codon_table, target_codons):
    """Compute relative codon usage within each amino acid group from absolute codon usage and exclude target codons"""
    
    #Add the sense codons to aa_groups
    aa_groups = {aa: [] for aa in codon_table.back_table if isinstance(aa, str)}
    for codon, aa in codon_table.forward_table.items():
        if codon not in target_codons:
            aa_groups[aa].append(codon)   
    
    # Add the stop codons to aa_groups
    aa_groups['Stop'] = []
    for codon in codon_table.stop_codons:
        if codon not in target_codons:
            aa_groups['Stop'].append(codon)
           
    #Compute the relative usage of all the non-target codons    
    relative_usage = {}
    for aa, codons in aa_groups.items():
        total = sum(codon_usage[codon] for codon in codons)
        for codon in codons:
            relative_usage[codon] = codon_usage[codon] / total
            
    return relative_usage

def get_CDS_locations(record):
    """Return a list of FeatureLocations for all CDS in a SeqRecord."""
    CDS_locations = []
    
    for feature in record.features:
        if feature.type == "CDS":

            if feature.location_operator == "join":
                for gene_part in feature.location.parts:
                    gene_part_location = gene_part
                    CDS_gene_part = Location(int(gene_part_location.start),int(gene_part_location.end))
                    CDS_locations.append(CDS_gene_part)
            else:
                element_location = Location(int(feature.location.start),int(feature.location.end))
                CDS_locations.append(element_location)            


    return CDS_locations

def compute_possible_codons(codon_table, aa, target_codons):
    """Compute possible codons for an amino acid, excluding target codons."""
    
    # Special case for stop codons
    if aa == 'Stop':
        possible_codons = codon_table.stop_codons
    else:
        # Create a list of all codons that code for the given amino acid
        possible_codons = [codon for codon, amino_acid in codon_table.forward_table.items() if amino_acid == aa]

    # Exclude the target codons from the list
    possible_codons = [codon for codon in possible_codons if codon not in target_codons]
    
    return possible_codons

def recode_codons(record, target_codons, restriction_sites, relative_usage):
    """Recode the codons in a DNA sequence"""
    
    recoded_seq = MutableSeq(record.seq)  # MutableSeq object
    codon_table = Bio.Data.CodonTable.standard_dna_table

    # Get CDS locations
    CDS_locations = get_CDS_locations(record)

    # Process each gene part
    for gene_part in CDS_locations:
        start, end = int(gene_part.start), int(gene_part.end)
        
        # Iterate over the gene part in codon-sized chunks
        for i in range(start, end, 3):
            codon = recoded_seq[i:i+3]
            
            # If the codon is a target, attempt to replace it
            if str(codon) in target_codons:
                if str(codon) in codon_table.stop_codons:
                    aa = 'Stop'
                else:
                    aa = codon_table.forward_table[str(codon)]

                possible_codons = compute_possible_codons(codon_table, aa, target_codons)
                # Prepare the probabilities for the synonymous codons
                weights = [relative_usage[codon] for codon in possible_codons]
                weights = [weight / sum(weights) for weight in weights] # normalize to sum to 1
                
                # Try synonymous codons based on codon usage bias
                for replacement in np.random.choice(possible_codons, size=60*len(possible_codons), p=weights):
                    if not would_create_restriction_site(recoded_seq, i, replacement, restriction_sites):
                        recoded_seq[i:i+3] = replacement
                        break


    #Convert MutableSeq back to SeqRecord with original metadata
    recoded_record = SeqRecord(Seq(recoded_seq), id=record.id, name=record.name,
                               description=record.description, annotations=record.annotations)
    recoded_record.features = record.features
    return recoded_record



def would_create_restriction_site(seq, i, replacement, restriction_sites):
    """Check if replacing a codon would create a restriction site"""
    
    # Replace codon and extract surrounding sequence
    temp_seq = seq[:i] + replacement + seq[i+3:]
    context = temp_seq[max(0, i-20):i+23]  # adjusted to 20-bp context
    
    # Search for restriction sites in the context
    for site in restriction_sites:
        if site in context:
            return True
    return False


def recode_restriction_sites(record, restriction_sites):
    
    sequence = str(record.seq)
    CDS_locations = get_CDS_locations(record)
    restriction_sites = list(restriction_sites)

    for enzyme_name in restriction_sites:
        # Get the recognition sequence for the enzyme
        enzyme = RestrictionBatch([enzyme_name])
        site = list(enzyme)[0].site

        # Find the locations of the restriction sites on both strands
        sites_fwd = [match.start() for match in re.finditer(site, sequence)]
        sites_rev = [match.start() for match in re.finditer(str(Seq(site).reverse_complement()), sequence)]
        site_locations = sites_fwd + sites_rev


        for loc in site_locations:
            # Check if the restriction site is in a CDS
            in_cds = any(location.start <= loc <= location.end for location in CDS_locations)
            site_end = loc + len(site)

            if not in_cds:
                # If the restriction site is not in a CDS, just optimize the site
                problem = DnaOptimizationProblem(
                    sequence=sequence,
                    constraints=[AvoidPattern(site, location=Location(loc, site_end))],
                    logger=None
                )

            elif sum(location.start <= loc <= location.end for location in CDS_locations) == 1:
                # If the restriction site is in exactly one CDS, optimize the site and enough surrounding bases to maintain the reading frame
                cds_location = [location for location in CDS_locations if location.start <= loc <= location.end][0]
                optimization_start = loc - (loc - cds_location.start) % 3
                optimization_end = site_end + (cds_location.end - site_end) % 3
                problem = DnaOptimizationProblem(
                    sequence=sequence,
                    constraints=[AvoidPattern(site, location=Location(optimization_start, optimization_end)), 
                                 EnforceTranslation(location=cds_location)],
                    logger=None
                )

            # Resolve constraints and replace the site in the original sequence with the optimized site
            try:
                problem.resolve_constraints()
                sequence = str(problem.to_record(with_sequence_edits=True).seq)
            except Exception as e:
                print(f"Could not optimize restriction site at location {loc}. Error: {e}")


    # Create a new SeqRecord with the optimized sequence and original metadata
    optimized_record = SeqRecord(Seq(sequence), id=record.id, name=record.name,
                             description=record.description, annotations=record.annotations)
    optimized_record.features = record.features

    return optimized_record


def write_to_genbank(record, filename):
    SeqIO.write([record], filename, "genbank")

def main():
    parser = argparse.ArgumentParser(description='Recode a phage genome by removing specified codons and restriction sites via synonymous mutations.')
    parser.add_argument('-i','--input', type=str, required=True, help='Input GenBank file name.')
    parser.add_argument('-o','--output', type=str, required=True, help='Output GenBank file name.')
    parser.add_argument('-s','--species', type=str, help='Species name for codon usage bias.')
    parser.add_argument('-c','--codons', type=str, nargs='+', help='Codons to recode.')
    parser.add_argument('-re','--re_sites', type=str, nargs='+', help='Restriction enzyme recognition sites to recode.')
    args = parser.parse_args()

    record = SeqIO.read(args.input, "genbank")
    codon_usage = get_codon_usage(args.species)
    restriction_sites = args.re_sites if args.re_sites else []


    if args.re_sites:
        record = recode_restriction_sites(record, restriction_sites)

    if args.codons:
        relative_usage = compute_relative_codon_usage(codon_usage, Bio.Data.CodonTable.standard_dna_table, args.codons)
        record = recode_codons(record, args.codons, restriction_sites, relative_usage)
        

    write_to_genbank(record, args.output)

if __name__ == "__main__":
    main()


