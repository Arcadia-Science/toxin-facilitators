#! /usr/bin/python3

from Bio import SeqIO
from Bio.SearchIO import parse

import subprocess
import argparse

"""
usage: run_all_pfam.py [-h] [--pfam_db PFAM_DB] [--evalue EVALUE] [--output_hmm OUTPUT_HMM] [--summary SUMMARY] fasta_file metadata_file

positional arguments:
  fasta_file            Protein FASTA file
  metadata_file         Metadata file of accessions to exclude

options:
  -h, --help            show this help message and exit
  --pfam_db PFAM_DB     Path to Pfam-A.hmm DB
  --evalue EVALUE       E-value cutoff for hmmsearch
  --output_hmm OUTPUT_HMM
                        Output HMMER file of Pfam hits
  --summary SUMMARY     Output txt file summarizing all Pfam hits excluding those in metadata file
"""

def run_hmmsearch(fasta_file, pfam_database, output_file, evalue=1e-50):
    """
    Run hmmsearch on input file against Pfam-A.hmm database
    """
    command = f"hmmsearch --tblout {output_file} -E {evalue} {pfam_database} {fasta_file}"
    subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def remove_hits(output_file, metadata_file, summary_file):
    """
    Parse the hmmsearch output file, remove hits that are in the metadata file
    """
    remove_ids = set()
    with open(metadata_file, "r") as metadata:
        for line in metadata:
            remove_ids.add(line.strip())

    # remove hits listed in metadata file from the output
    filtered_hits = {}
    for qresult in parse(output_file, "hmmer3-tab"):
        query_accession = qresult.accession
        if query_accession not in remove_ids:
            for hit in qresult.hits:
                filtered_hits[hit.id] = {
                    "query_accession": query_accession,
                    "query_name": qresult.id
                }

    # write summary of filtered hits to file
    with open(summary_file, "w") as output:
        output.write("hit_id\tpfam_query_accession\tpfam_name\n")
        for hit_id in filtered_hits:
            hit=filtered_hits[hit_id]
            output.write(f"{hit_id}\t{hit['query_accession']}\t{hit['query_name']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_file", help="Protein FASTA file")
    parser.add_argument("metadata_file", help="Metadata file of accessions to exclude")
    parser.add_argument("--pfam_db", help="Path to Pfam-A.hmm DB")
    parser.add_argument("--evalue", type=float, default=1e-50, help="E-value cutoff for hmmsearch")
    parser.add_argument("--output_hmm", default="all_pfam_hits.out", help="Output HMMER file of Pfam hits")
    parser.add_argument("--summary", default="all_additional_pfam_annotations.txt", help="Output txt file summarizing all Pfam hits excluding those in metadata file")
    args = parser.parse_args()

    run_hmmsearch(args.fasta_file, args.pfam_db, args.output_hmm, args.evalue)
    remove_hits(args.output_hmm, args.metadata_file, args.summary)
