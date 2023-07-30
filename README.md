# Bacteriophage Genome Recoding Script

## Description

This Python script is designed to recode bacteriophage genomes by replacing specific codons and restriction enzyme recognition sites with synonymous codons. Synonymous codons are chosen based on their optimized usage for a specified host organism.

The script uses the Biopython library to read and manipulate the GenBank file, compute the codon usage of the host species, and recode the bacteriophage genome. After recoding, the script checks for any changes in protein translations, identifies restriction sites and codons in the new genome that were not recoded, and writes the recoded genome to a new GenBank file. Finally, it generates a report detailing the results of the recoding process.

## Dependencies

This script requires the following Python libraries:

- Biopython
- DNAChisel
- Numpy
- Pandas
- Argparse
- Regular Expressions (re)
- tqdm
- datetime
- os

## How to Use

1. Install the necessary Python packages.

Either:

'pip install biopython dnachisel numpy pandas tqdm'

Or: 'conda install -c conda-forge biopython dnachisel numpy pandas tqdm'

2. Clone the repository to your local machine.

'git clone https://github.com/petevoor/genome_recoder.git'

3. Navigate to the directory containing the script.

'cd directory_name'

4. Run the script with necessary arguments.

'python recode_genome.py -i input_file.gb -o output_file.gb -s species_name -c codon1 codon2 -re enzyme1 enzyme2'

### Arguments

- `-i` or `--input`: Input GenBank file name. Specify as path name if input file is not in the same directory as recode_genome.py.
- `-o` or `--output`: Output GenBank file name. Specify as path name to specify output directory.
- `-s` or `--species`: Species name for codon usage bias.
- `-c` or `--codons`: Space-separated list of codons to recode. E.g., TCA TCG.
- `-re` or `--re_sites`: Space-separated list of restriction enzyme recognition sites to recode. E.g., BsaI BsmbI.

## Limitations

- This script assumes that the input genome is correctly annotated and that the GenBank file contains all the necessary information for the recoding process.
- Currently, this script is unable to recode codons or restriction sites that fall in overlapping, out of frame coding sequences.
- The codon usage bias for the specified species must be available on the Codon Usage Database.
- The efficiency of recoding might vary based on the chosen host organism and the codon bias.
- Recoding might introduce unintended changes to the proteome, and it is always recommended to check the recoding report and GenBank file after running the script.