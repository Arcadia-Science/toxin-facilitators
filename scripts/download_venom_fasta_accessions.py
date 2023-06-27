#! /usr/bin/python3

import argparse
import subprocess
import os
import time

"""
usage: download_venom_fasta_accessions.py [-h] [-l LOG_FILE] [-f FORMAT] [-d DELAY] metadata_file

Download fasta sequences from NCBI for each accession in the metadata file

positional arguments:
  metadata_file         Path to metadata CSV file

optional arguments:
  -h, --help            show this help message and exit
  -l LOG_FILE, --log_file LOG_FILE
  -f FORMAT, --format FORMAT
                        Format for efetch command (default: fasta)
  -d DELAY, --delay DELAY
                        Delay between each download in seconds (default: 180)
"""

def create_accession_dict(metadata_file):
    accession_dict = {}
    # read in header lines
    with open(metadata_file, "r", encoding="utf-8-sig") as f:
        header = f.readline().strip().split(',')
        accession_id = header.index('accession')
        species_id = header.index('species_name')
        contigs_id = header.index('contigs')
    # create dictionaries
        for line in f:
            fields = line.strip().split(',')
            accession = fields[accession_id]
            species_name = fields[species_id].replace(' ', '_')
            contigs = int(fields[contigs_id])
            accession_dict[accession] = {'species_name': species_name, 'contigs': contigs}
    return accession_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Download fasta sequences from NCBI for each accession in the metadata file")
    parser.add_argument('metadata_file', type=str, help="Path to metadata CSV file")
    parser.add_argument('-l', '--log_file', type=str, default="download_status.log")
    parser.add_argument('-f', '--format', type=str, default="fasta", help="Format for efetch command (default: fasta)")
    parser.add_argument('-d', '--delay', type=int, default=180, help="Delay between each download in seconds (default: 180)")
    args = parser.parse_args()

    accession_dict = create_accession_dict(args.metadata_file)

    with open(args.log_file, 'a') as log:
        for accession, info in accession_dict.items():
            species_name = info['species_name']
            # create list of accession codes to all pass to efetch since it can take comma-separated accessions in one ping
            accession_codes = []
            total_contigs = info['contigs']
            for contig_num in range(1, total_contigs + 1):
                accession_code = f"{accession}{contig_num:06}"
                accession_codes.append(accession_code)

            # pass to efetch and write to species outfile
            outfile = f"{species_name}_{accession}.fasta"
            efetch_cmd = ["efetch", "-db", "nuccore", "-id", ",".join(accession_codes), "-format", args.format]
            with open(outfile, "a") as f:
                efetch_process = subprocess.Popen(efetch_cmd, stdout=f)


            # delay between each species to pull down from
            time.sleep(args.delay)

            # check if output file is empty and write to a log file the status of each species/accession
            if os.path.getsize(outfile) == 0:
                os.remove(outfile)
                log.write(f"Failed to download fasta accessions {accession} for {species_name}!\n")
            else:
                log.write(f"Downloaded fasta accessions {accession} for {species_name} to {outfile} successfully!\n")
