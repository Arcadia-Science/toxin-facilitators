#! /usr/bin/python3

from Bio import SeqIO

import csv
import argparse

def read_metadata(metadata_file):
    metadata = {}
    with open(metadata_file, "r") as input:
        header = input.readline().strip().split('\t')
        entry_id = header.index('Entry')
        protein_name = header.index('Protein names')
        for line in input:
            fields = line.strip().split('\t')
            accession = fields[entry_id]
            description = fields[protein_name]
            metadata[accession] = description
    return metadata

def get_protein_info(fasta_file, metadata_tsv, results_tsv, output_file):
    # first read in metadata tsv to run function and store in dictionary
    metadata = read_metadata(metadata_tsv)
    # open the output results file
    with open(output_file, "w") as outfile:
        writer = csv.writer(outfile, delimiter = "\t")
        writer.writerow(["protein_id", "length", "cysteine_content", "accession_hit", "sequence_identity", "alignment_length", "evalue", "bit_score", "accession_hit_description"])
    # open results tsv from mmseqs
        with open(results_tsv, "r") as input:
            reader = csv.reader(input, delimiter = "\t")
            # read in fasta file and store all record in dictionary indexed by the ID to look up the protein_id
            records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

        # loop over rows for every protein id to get additional needed columns
            for row in reader:
                accession = row[0]
                protein_id = row[1]
                sequence_identity = row[2]
                alignment_length = row[3]
                evalue = row[10]
                bit_score = row[11]
                # look up protein ID in the records dictionary
                record = records.get(protein_id)
                if record is None:
                    print(f"Protein ID {protein_id} not found in FASTA file")
                    continue
                sequence = str(record.seq)
                length = len(sequence)
                cysteine_content = sequence.count("C") / length
                description = metadata.get(accession, "")
                writer.writerow([protein_id, length, cysteine_content, accession, sequence_identity, alignment_length, evalue, bit_score, description])

def parse_arguments():
    parser = argparse.ArgumentParser(description="Get protein information from a TSV file and a FASTA file")
    parser.add_argument("results_tsv", help="path to TSV file containing protein IDs")
    parser.add_argument("fasta_file", help="path to FASTA file containing protein sequences")
    parser.add_argument("metadata_tsv", help="path to TSV file containing metadata for accessions")
    parser.add_argument("output_file", help="path to output file for protein information")
    return parser.parse_args()

args = parse_arguments()

get_protein_info(args.fasta_file, args.metadata_tsv, args.results_tsv, args.output_file)
