#! /usr/bin/python3

from Bio import SeqIO
from Bio.SearchIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

import argparse
import subprocess

"""
usage: run_parse_hmms.py [-h] [--evalue EVALUE] [--output_file OUTPUT_FILE] hmm_file fasta_file

positional arguments:
  hmm_file              HMM profiles file
  fasta_file            Protein FASTA file

options:
  -h, --help            show this help message and exit
  --evalue EVALUE       E-value cutoff for hmmsearch
  --output_file OUTPUT_FILE
                        Output FASTA file for hits
"""

def run_hmmsearch(hmm_file, fasta_file, evalue=1e-50):
    """
    Run hmmsearch on given input FASTA file and HMM profile
    """
    command = f"hmmsearch --tblout hits.txt -E {evalue} {hmm_file} {fasta_file}"
    subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def parse_hits_output(hmmer_output):
    """
    Parse the hits file and return list of locus tag hits
    """
    hits_dict = {}
    for qresult in parse(hmmer_output, "hmmer3-tab"):
        query_name  = qresult.id
        query_accession = qresult.accession
        for hit in qresult.hits:
            hit_id = hit.id
            evalue = hit.hsps[0].evalue
            hits_dict[hit_id] = (query_accession, query_name, evalue)
    return hits_dict

    return hits_dict

def write_fasta_file(hits, input_fasta, hits_fasta):
    """
    Write hits to a new fasta file
    """
    records = list(SeqIO.parse(input_fasta, "fasta"))
    with open(hits_fasta, "a") as f:
        for hit in hits:
            for record in records:
                if record.id == hit:
                    f.write(f">{hit}\n")
                    f.write(f"{record.seq}\n")

def summarize_hits(hits_dict, hits_fasta, summary_file):
    """
    Summarize protein length, cysteine content, the hit description
    """
    with open(summary_file, "w") as f:
        f.write("hit_id\tlength\tcysteine_content\te_value\tquery_accession\tquery_name\n")
        for hit_id, (query_accession, query_name, evalue) in hits_dict.items():
            for record in SeqIO.parse(hits_fasta, "fasta"):
                if record.id == hit_id:
                    length = len(record.seq)
                    cysteine_count = record.seq.count("C") + record.seq.count("c")
                    cysteine_content = cysteine_count / length if length > 0 else 0
                    f.write(f"{hit_id}\t{length}\t{cysteine_content:.2f}\t{evalue}\t{query_accession}\t{query_name}\n")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("hmm_file", help="HMM profile file")
    parser.add_argument("fasta_file", help="Protein FASTA file")
    parser.add_argument("--evalue", type=float, default=1e-50, help="E-value cutoff for hmmsearch")
    parser.add_argument("--output_file", default="hits.fasta", help="Output FASTA file for hits")
    args = parser.parse_args()

    # Run hmmsearch and parse the results
    run_hmmsearch(args.hmm_file, args.fasta_file, args.evalue)
    hits_dict = parse_hits_output("hits.txt")

    # Write the hits to a single FASTA file
    write_fasta_file(list(hits_dict.keys()), args.fasta_file, args.output_file)

    # summarize hits directly from the fasta file and write to a summary file
    summarize_hits(hits_dict, args.output_file, "summary.txt")
