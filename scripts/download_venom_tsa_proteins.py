#! /usr/bin/python3

import argparse
import subprocess
import os
import time

"""
usage: download_venom_tsa_proteins.py [-h] [-l LOG_FILE] [-f FORMAT] [-d DELAY] metadata_file

Download protein sequences from NCBI for each accession in the metadata file

positional arguments:
  metadata_file         Path to metadata CSV file

optional arguments:
  -h, --help            show this help message and exit
  -l LOG_FILE, --log_file LOG_FILE
  -f FORMAT, --format FORMAT
                        Format for efetch command (default: fasta_cds_aa)
  -d DELAY, --delay DELAY
                        Delay between each download in seconds (default: 180)
"""

def create_accession_to_species_dict(metadata_file):
    accession_to_species = {}
    # read in header lines
    with open(metadata_file, 'r', encoding='utf-8-sig') as f:
        header = f.readline().strip().split(',')
        accession_id = header.index('accession')
        species_id = header.index('species_name')
    # read in remaining lines and create dictionary
        for line in f:
            fields = line.strip().split(',')
            accession = fields[accession_id]
            species_name = fields[species_id].replace(' ','_')
            accession_to_species[accession] = species_name
    return accession_to_species

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Download protein sequences from NCBI for each accession in the metadata file")
    parser.add_argument('metadata_file', type=str, help="Path to metadata CSV file")
    parser.add_argument('-l', '--log_file', type=str, default="download_status.log")
    parser.add_argument('-f', '--format', type=str, default="fasta_cds_aa", help="Format for efetch command (default: fasta_cds_aa)")
    parser.add_argument('-d', '--delay', type=int, default=180, help="Delay between each download in seconds (default: 180)")
    args = parser.parse_args()

    accession_to_species = create_accession_to_species_dict(args.metadata_file)

    with open(args.log_file, 'a') as log:
        for accession, species_name in accession_to_species.items():

            # call esearch and efetch separately and write to outfile
            # this seemed to be one of the only ways that this would behave correctly and download successfully
            outfile = f"{species_name}_{accession}.fasta"
            esearch_cmd = ["esearch", "-db", "nuccore", "-query", accession]
            efetch_cmd = ["efetch", "-format", args.format]

            with open(outfile, 'w') as f:
                esearch_process = subprocess.Popen(esearch_cmd, stdout=subprocess.PIPE)
                efetch_process = subprocess.Popen(efetch_cmd, stdin=esearch_process.stdout, stdout=f)
                esearch_process.stdout.close()
                efetch_process.communicate()

            # delay between each download
            time.sleep(args.delay)

            # check if output file is empty and write to a log file the status of each species/accession
            if os.path.getsize(outfile) == 0:
                os.remove(outfile)
                log.write(f"Failed to download protein accession {accession} for {species_name}!\n")
            else:
                log.write(f"Downloaded protein accession {accession} for {species_name} to {outfile} successfully!\n")
