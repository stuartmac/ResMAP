"""
Counts the number of occurrences of each sequence in a library in a BAM file.

Authors: 
    Stuart MacGowan <smacgowan@dundee.ac.uk>
    Richard Wall <Richard.Wall@lshtm.ac.uk>

Date:
    2023-12-16
"""

import argparse
import csv

import pysam

from tqdm import tqdm


# Parse command line arguments
parser = argparse.ArgumentParser(description='Counts sequences in library')
parser.add_argument('-i', nargs='?', dest='bam', help='Input BAM file')
parser.add_argument('-o', nargs='?', dest='output', help='Output TSV file')
parser.add_argument('--library', nargs='?', dest='library', help='Library csv file')
args = parser.parse_args()

# Validate input arguments
if not args.bam or not args.library:
    parser.error("Both BAM file and library CSV file are required.")
    exit(1)

# Load library
library = {}
with open(args.library, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    for row in reader:
        library[row[0]] = row[1]

# Count library sequences
counts = {}
with pysam.AlignmentFile(args.bam, "rb") as bam:
    # Wrap bam.fetch() with tqdm for the progress bar
    for read in tqdm(bam.fetch("Pf3D7_13_v3", 240, 310), desc="Processing BAM file", unit=" reads"):
        # if read.is_unmapped:
        #     continue
        if read.reference_start < 160:
            continue
        for label, seq in library.items():
            if seq in read.query_sequence:
                if label not in counts:
                    counts[label] = 0
                counts[label] += 1
                # Break out of the loop after the first match
                break

# Write counts to tsv
output_file = args.output if args.output else "counts.tsv"
try:
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(['Label', 'Sequence', 'Count'])
        for label, count in counts.items():
            writer.writerow([label, library[label], count])
    print(f"Counts written to {output_file}")
except IOError as e:
    print(f"Error writing to file {output_file}: {e}")
