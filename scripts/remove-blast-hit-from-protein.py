#! /usr/bin/python3

from Bio import SeqIO
import argparse


# arguments
parser = argparse.ArgumentParser(description='Remove parts of protein sequences based on start and end positions.')
parser.add_argument('tsv_file', help='TSV file from BLAST results containing protein locus tag, start position, and end position')
parser.add_argument('input_fasta_file', help='input FASTA file')
parser.add_argument('output_fasta_file', help='output FASTA file')

args = parser.parse_args()

# read in a TSV file into a dictionary with locus tag as key, and list of start/end positions
protein_ranges = {}
with open(args.tsv_file, "r") as tsv_file:
    for i, line in enumerate(tsv_file):
        locus_tag, start, end = line.strip().split('\t')
        if locus_tag not in protein_ranges:
            protein_ranges[locus_tag] = []
        protein_ranges[locus_tag].append((int(start), int(end)))

# iterate through input fasta and write the modified sequences
seen_sequences = {} # hold the sequences to only add to the fasta file if that exact sequence hasn't been added before, removes duplicates from the exact same species range from the duplicate BLAST hit problem
with open(args.output_fasta_file, "w") as output_fasta_file:
    with open(args.input_fasta_file, "r") as input_fasta_file:
        for record in SeqIO.parse(input_fasta_file, "fasta"):
            if record.id in protein_ranges:
                for i, (start, end) in enumerate(protein_ranges[record.id]):
                    modified_sequence = record.seq[:start-1] + record.seq[end:]
                    if len(modified_sequence) > 1: # only remove empty records, keep things even if short
                        sequence_str = str(modified_sequence)
                        if sequence_str not in seen_sequences:
                            modified_record = record[:]
                            modified_record.seq = modified_sequence
                            if i == 0:
                                modified_record.id = record.id
                            else:
                                modified_record.id = f"{record.id}_rep{i+1}" # names replicate records that have slightly different sequence outputs
                            SeqIO.write(modified_record, output_fasta_file, "fasta")
                            seen_sequences[sequence_str] = modified_record.id
                        else:
                            if i == 0:
                                seen_sequences[sequence_str] = record.id
                            else:
                                seen_sequences[sequence_str] = f"{record.id}_rep{i+1}"
