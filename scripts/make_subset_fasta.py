#! /usr/bin/python3


from Bio import SeqIO
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Subset a FASTA file based on protein locus tags.')
parser.add_argument('input_locus_tags', help='File with list of protein locus tags')
parser.add_argument('input_fasta', help='input FASTA file')
parser.add_argument('output_fasta', help='output FASTA file')

args = parser.parse_args()

# read in protein locus tags
with open(args.input_locus_tags, "r") as locus_tags:
    locus_tag_set = set(line.strip() for line in locus_tags)

# iterate through fasta records, for each matching locus tags write to a new fasta file
with open(args.output_fasta, "w") as output_fasta_file:
    with open(args.input_fasta, "r") as input_fasta_file:
        for record in SeqIO.parse(input_fasta_file, "fasta"):
            if record.id in locus_tag_set:
                SeqIO.write(record, output_fasta_file, "fasta")
