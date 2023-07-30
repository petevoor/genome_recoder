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
from tqdm import tqdm
from datetime import date
import os



# Load the CUTG csv file into a DataFrame
df = pd.read_csv('dbs/dna_codon_usage.csv')

recoded_sites_codons = []

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

def find_overlapping_regions(record):
    """Find regions in the genome where two or more CDS overlap."""
    CDS_locations = get_CDS_locations(record)

    overlaps = set()
    for i in range(len(CDS_locations)):
        for j in range(i + 1, len(CDS_locations)):
            if (CDS_locations[i].end > CDS_locations[j].start and
                CDS_locations[i].start < CDS_locations[j].end):
                # The overlapping region is from the start of the next CDS to the end of the current CDS
                overlap_start = max(CDS_locations[i].start, CDS_locations[j].start)
                overlap_end = min(CDS_locations[i].end, CDS_locations[j].end)
                # If the overlapping regions are not in the same reading frame
                if overlap_start % 3 != CDS_locations[i].start % 3:
                    overlaps.add((overlap_start, overlap_end))

    return overlaps


def would_create_restriction_site(seq, i, replacement, restriction_sites):
    """Check if replacing a codon would create a restriction site"""
    
    # Replace codon and extract surrounding sequence
    temp_seq = seq[:i] + replacement + seq[i+3:]
    context = temp_seq[max(0, i-20):i+23]  # adjusted to 20-bp context

    # Search for restriction sites in the context
    for enzyme_name in restriction_sites:
        # Get the recognition sequence for the enzyme
        enzyme = RestrictionBatch([enzyme_name])
        site = list(enzyme)[0].site
        if site in context or str(Seq(site).reverse_complement()) in context:
            return True
    return False


def recode_codons(record, target_codons, restriction_sites, relative_usage):
    """Recode the codons in a DNA sequence"""
    
    recoded_seq = MutableSeq(record.seq)  # MutableSeq object
    codon_table = Bio.Data.CodonTable.standard_dna_table
    
    # Find overlapping regions
    overlaps = find_overlapping_regions(record)

    # Get CDS locations
    CDS_locations = get_CDS_locations(record)

    # Process each gene part
    for gene_part in tqdm(CDS_locations, desc='Recoding codons', unit='CDS'):
        start, end = int(gene_part.start), int(gene_part.end)
        
        # Iterate over the gene part in codon-sized chunks
        for i in range(start, end, 3):
            
            # If the codon is in an overlapping region, skip it
            if any(overlap_start <= i < overlap_end for overlap_start, overlap_end in overlaps):
                continue
            
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
                        recoded_sites_codons.append((i, str(codon), replacement))
                        break


    #Convert MutableSeq back to SeqRecord with original metadata
    recoded_record = SeqRecord(Seq(recoded_seq), id=record.id, name=record.name,
                               description=record.description, annotations=record.annotations)
    recoded_record.features = record.features
    return recoded_record


def recode_restriction_sites(record, restriction_sites):
    
    sequence = str(record.seq)
    CDS_locations = get_CDS_locations(record)
    restriction_sites = list(restriction_sites)
    
    # Find overlapping regions
    overlaps = find_overlapping_regions(record)

    for enzyme_name in tqdm(restriction_sites, desc='Removing restriction sites', unit='site'):
        # Get the recognition sequence for the enzyme
        enzyme = RestrictionBatch([enzyme_name])
        site = list(enzyme)[0].site

        # Find the locations of the restriction sites on both strands
        sites_fwd = [match.start() for match in re.finditer(site, sequence)]
        sites_rev = [match.start() for match in re.finditer(str(Seq(site).reverse_complement()), sequence)]
        site_locations = sites_fwd + sites_rev


        for loc in site_locations:
            # If the restriction site is in an overlapping region, skip it
            if any(overlap_start <= loc < overlap_end for overlap_start, overlap_end in overlaps):

                continue
            
            # Check if the restriction site is in a CDS
            in_cds = any(location.start <= loc <= location.end for location in CDS_locations)
            site_end = loc + len(site)

            if not in_cds:
                continue

            elif any(location.start <= loc <= location.end for location in CDS_locations):

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
                recoded_sites_codons.append((loc, site, sequence[loc:loc+len(site)]))
            except Exception as e:
                print(f"Could not optimize restriction site at location {loc}. Error: {e}")



    # Create a new SeqRecord with the optimized sequence and original metadata
    optimized_record = SeqRecord(Seq(sequence), id=record.id, name=record.name,
                             description=record.description, annotations=record.annotations)
    optimized_record.features = record.features

    return optimized_record


def write_to_genbank(record, output_path, codons, restriction_sites):
    record.id = "recoded_" + record.id
    
    description = record.description
    if codons:
        description += f" (this organism has been recoded to remove {' '.join(codons)} codons"
        if restriction_sites:
            description += " and"
    if restriction_sites:
        description += f" {' '.join(restriction_sites)} restriction enzyme recognition sites"
    description += f" from the genome of {record.annotations.get('source', '')})"
    
    record.description = description
    SeqIO.write(record, output_path, "genbank")
    

def check_translations(ref_record, new_record):
    """Check that the translations of the CDS in the old and new records are identical."""
    # Initialize a list to store results
    discrepancies = {}
    reference = SeqIO.read(ref_record, "genbank")

    # Iterate over the features in the old and new records
    for old_feature, new_feature in zip(reference.features, new_record.features):
        # Only consider features that are CDS
        if old_feature.type == 'CDS' and new_feature.type == 'CDS':
            old_seq = old_feature.extract(reference.seq)
            new_seq = new_feature.extract(new_record.seq)
            
            # Translate the sequences
            old_protein = old_seq.translate(to_stop=True)
            new_protein = new_seq.translate(to_stop=True)
            
            # Check if the translations are the same
            if old_protein != new_protein:
                protein_name = new_feature.qualifiers.get("note")[0]
                discrepancies[protein_name] = {}
                for idx, (ref_aa, query_aa) in enumerate(zip(old_protein, new_protein)):
                    if ref_aa != query_aa:
                        discrepancies[protein_name][idx + 1] = [ref_aa, query_aa]

                         
    return discrepancies

def find_unrecoded_restriction_sites(new_record, restriction_sites):
    new_sequence = str(new_record.seq)
    unrecoded_restriction_sites = []
    
    for enzyme_name in restriction_sites:
        enzyme = RestrictionBatch([enzyme_name])
        site = list(enzyme)[0].site
        new_site_locations_fwd = [match.start() for match in re.finditer(site, new_sequence)]
        new_site_locations_rev = [match.start() for match in re.finditer(str(Seq(site).reverse_complement()), new_sequence)]
        new_site_locations = new_site_locations_fwd + new_site_locations_rev
        
        for loc in new_site_locations:
            unrecoded_restriction_sites.append((loc, site))
    
    return unrecoded_restriction_sites


def find_unrecoded_codons(ref_record, new_record, target_codons):
    reference = SeqIO.read(ref_record, "genbank")
    original_sequence = str(reference.seq)
    new_sequence = str(new_record.seq)
    unrecoded_codons = []
    
    # Get CDS locations
    CDS_locations = get_CDS_locations(reference)

    # Process each gene part
    for gene_part in CDS_locations:
        start, end = int(gene_part.start), int(gene_part.end)
        
        # Iterate over the gene part in codon-sized chunks
        for i in range(start, end, 3):
            original_codon = original_sequence[i:i+3]
            new_codon = new_sequence[i:i+3]
            
            # If the codon was a target and it is the same in the new sequence, it was not recoded
            if original_codon in target_codons and original_codon == new_codon:
                unrecoded_codons.append((i, original_codon))

    return unrecoded_codons


def generate_report(proteome_QC, unrecoded_restriction_sites, unrecoded_codons, recoded_sites_codons, filename):
    
    # Today's date
    today = date.today()
    
    # Removing extension from filename
    filename_without_ext = os.path.splitext(filename)[0]
    
    # File name for the report
    report_filename = str(today)+"_"+str(filename_without_ext)+'_Recoding_Report.csv'

    # Initialize the report string with title, date and overview
    report = f"{filename_without_ext} Recoding Report\nDate: {today}\n\n"
    
    # Number of each codon and of each restriction site that were successfully recoded
    num_recoded_codons = sum([1 for index, old_seq, new_seq in recoded_sites_codons if len(old_seq) == 3 and old_seq != new_seq])
    num_recoded_sites = sum([1 for index, old_seq, new_seq in recoded_sites_codons if len(old_seq) != 3 and old_seq != new_seq])

    report += f"Number of codons recoded: {num_recoded_codons}\n"
    report += f"Number of restriction sites recoded: {num_recoded_sites}\n"
    
    # Overview of whether there are any mutations to the proteome and, if so, how many
    num_mutations = len(proteome_QC)
    report += f"Number of mutations to the proteome: {num_mutations}\n"

    # Open the report file in write mode
    with open(report_filename, 'w') as report_file:
        # Write the initialized report string
        report_file.write(report)
        
        # If there are mutations to the proteome, write the proteome_QC table to the report
        if num_mutations > 0:
            report_file.write("\nProteome mutations:\n")
            proteome_QC.to_csv(report_file, sep=',',index=False)
        
        # Write the table of restriction sites that couldn't be recoded
        report_file.write(f"\nUnrecoded restriction sites ({len(unrecoded_restriction_sites)}):\n")
        unrecoded_restriction_sites_df = pd.DataFrame(unrecoded_restriction_sites, columns=['Index', 'Restriction Site'])
        unrecoded_restriction_sites_df.to_csv(report_file, sep=',',index=False)

        # Write the table of codons that couldn't be recoded
        report_file.write(f"\nUnrecoded codons ({len(unrecoded_codons)}):\n")
        unrecoded_codons_df = pd.DataFrame(unrecoded_codons, columns=['Index', 'Codon'])
        unrecoded_codons_df.to_csv(report_file, sep=',',index=False)

        # Write the table of all recoded restriction sites and codons
        report_file.write("\nRecoded sites and codons:\n")
        recoded_sites_codons_df = pd.DataFrame(recoded_sites_codons, columns=['Index', 'Old Sequence', 'New Sequence'])
        recoded_sites_codons_df.to_csv(report_file, sep=',',index=False)
        
    print(f"Report written to {report_filename}")


def main():
    parser = argparse.ArgumentParser(description='Recode a phage genome by removing specified codons and restriction sites via synonymous mutations.')
    parser.add_argument('-i','--input', type=str, required=True, help='Input GenBank file name. Specify as path name if input file is not in the same directory as recode_genome.py')
    parser.add_argument('-o','--output', type=str, required=True, help='Output GenBank file name. Specify as path name to specify output directory')
    parser.add_argument('-s','--species', type=str, help='Species name for codon usage bias.')
    parser.add_argument('-c','--codons', type=str, nargs='+', help='Space separated list of codons to recode. e.g. TCA TCG')
    parser.add_argument('-re','--re_sites', type=str, nargs='+', help='Space separated list of restriction enzyme recognition sites to recode. e.g. BsaI BsmbI')
    args = parser.parse_args()

    record = SeqIO.read(args.input, "genbank")
    codon_usage = get_codon_usage(args.species)
    restriction_sites = args.re_sites if args.re_sites else []


    if args.re_sites and args.codons:
        record = recode_restriction_sites(record, restriction_sites)
        relative_usage = compute_relative_codon_usage(codon_usage, Bio.Data.CodonTable.standard_dna_table, args.codons)
        record = recode_codons(record, args.codons, restriction_sites, relative_usage)
        unrecoded_restriction_sites = find_unrecoded_restriction_sites(record, restriction_sites)
        unrecoded_codons = find_unrecoded_codons(args.input, record, args.codons)
    
    elif args.re_sites:
        record = recode_restriction_sites(record, restriction_sites)
        unrecoded_restriction_sites = find_unrecoded_restriction_sites(record, restriction_sites)
        unrecoded_codons = []
    
    elif args.codons:
        relative_usage = compute_relative_codon_usage(codon_usage, Bio.Data.CodonTable.standard_dna_table, args.codons)
        record = recode_codons(record, args.codons, restriction_sites, relative_usage)
        unrecoded_restriction_sites = []
        unrecoded_codons = find_unrecoded_codons(args.input, record, args.codons)

        
    proteome_QC = check_translations(args.input, record)
    
    write_to_genbank(record, args.output, args.codons, restriction_sites)
    

    if proteome_QC:
        print("Proteome mutations found in the following locations:")
        for protein_name, mutation_info in proteome_QC.items():
            print(f"Protein: {protein_name}")
            for location, (old_aa, new_aa) in mutation_info.items():
                print(f"Location: {location}, Old Amino Acid: {old_aa}, New Amino Acid: {new_aa}")
    else:
        print("No mutations found in the proteome.")
        
    # Generate report after all processes are done
    generate_report(proteome_QC, unrecoded_restriction_sites, unrecoded_codons, recoded_sites_codons, args.output)
    print(len(unrecoded_codons))
    print(len(unrecoded_restriction_sites))
    print(unrecoded_restriction_sites)

if __name__ == "__main__":
    main()


