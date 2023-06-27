# Discovery of Novel Toxin Facilitators
The strike team is interested in finding non-toxic elements that may facilitate toxins in venoms or poisons reach their final destination. This could encompass proteins, peptides, small molecules, lipids, etc. Mining proteomes from diverse venoms might lead to promising non-toxic elements that can aid in drug delivery.

This README first summarizes the computational approach we are taking for mining venom proteomes for facilitator proteins that might aid toxins to persist in the host environment and reach their target destination.

## Background
- We hypothesize that there are conserved facilitator proteins that are non-toxic but facilitate toxins in venoms and poisons to reach their final destination
- The first two ideas under the umbrella of this hypothesis are that toxins and potentially facilitator proteins are 1) Cysteine rich (since toxins have a lot of disulfide bonds) and 2) Potentially are fused to the toxin itself.
- There are several publicly available transcriptomes of venom from animals, insects, etc. and can be used for _de novo_ transcriptome assembly and annotation to obtain proteins or directly retrieve the proteins already annotated from these datasets and perform downstream analyses.

## Planned Approach
1. Download proteins from the TSA for venom gland accessions that have proteins available. Otherwise download transcribed RNA sequences and predict ORFs/proteins, currently using `transdecoder`.
2. Search venom proteins against HMM profiles of known toxins from the Pfam domain database to filter for proteins with domains that look like toxins and annotate toxin+ looking candidate proteins with all Pfam HMMs
3. Perform clustering of representative sequences with `mmseqs` for curated toxin proteins in Uniprot [https://www.uniprot.org/help/Toxins](https://www.uniprot.org/help/Toxins)
4. Explore the resulting candidates taking into account cysteine content, protein length, domains that look like toxins. If the protein has an additional domain other than the original hit to the toxin database, explore the sequence/structure space of that domain.
5. Match up these results with protein folding analyses the annotation team is working on.

## Implementation
1. The first part of setting up this project involves downloading proteins/transcribed RNA/transcriptomes from venom gland entries in the TSA, some babysitting involved
2. The second part takes input proteins and a concatenated file of target search HMMs to run `hmmsearch` and filter based on TC or e-value. Also use this to propagate summary statistics such as cysteine content and protein length. Then, those hit proteins and annotates against the entire Pfam HMM database to get additional annotations besides the toxins themselves. This would help us find potential "fusion" events if there are multiple annotations/domains in a single protein
3. An additional but related analysis was added after #2 didn't quite work out since the Pfam toxin HMM domains are pretty general. Instead we are taking the full-length protein sequences of a curated set of toxins from Uniprot [https://www.uniprot.org/help/Toxins](https://www.uniprot.org/help/Toxins) and clustering into representative sequences with `mmseqs`, then searching against the venom proteins with `mmseqs easy-search`
4. After the clustering-based approach against the Uniprot toxin database, we discovered there were outlier proteins that were longer than the reference protein, and in some cases this pointed to proteins with additional domains such as mucins. We want to explore outlier proteins with additional domains in a more high-throughput way, since the manual way involved grabbing the non-toxin part of the protein by editing the FASTA file and BLASTing against the NCBI-nr database. Outliers in each cluster were defined as a mild outlier - anything Q3+1.5*IQ and pulled out for further analysis.

Listed below are instructions for setting up individual conda environments for the different pieces of software used to accomplish the above tasks. You can either do that or use the `environment.yml` file to create a conda environment for all software used for this project with:

```
conda env create -n toxin_facilitators -f environment.yml
```

### tl;dr of below
1. Download proteins or transcribed rna accessions with the `download_venom_tsa_proteins.py` or `download_venom_fasta_accessions.py` scripts.
2. Prepare the venom HMM database from the full Pfam database with the `metadata/toxin_pfam_accessions.txt` file, search the venom proteins against toxin Pfam HMMs with `run_parse_hmms.py`, then add additional Pfam annotations with `run_all_pfam.py`. These results are parsed with `toxin_facilitator_report.md`
3. Use an `mmseqs` clustering approach to first obtain representative sequences from the Uniprot curated set of toxin proteins, use those to run `mmseqs search` against the reference venom proteins, and get sets of toxin "clusters" to summarize
4. From a list of outlier proteins that are curated from the clusters, align against a DIAMOND database of the venom proteins and extract out the toxin part of those proteins to only have the accessory part leftover. Then either annotate with Pfam HMMs or BLAST to the NCBI-nr database.

 ## Downloading and preparing data

 The first part of this project that will be used for all downstream steps including `hmmsearch` against Pfam toxin domains and `mmseq easy-search` clustering against full-length Uniprot toxin proteins involves downloading protein accessions or transcribed RNA files from venom gland Transcriptome Shotgun Assemblies (TSA) datasets that are publicly available.

### Download proteins from the Transcriptome Shotgun Assemblies (TSAs) for venom glands
First install the Entrez Direct toolkit from NCBI following these instructions: [https://www.ncbi.nlm.nih.gov/books/NBK179288/](https://www.ncbi.nlm.nih.gov/books/NBK179288/). If you need help finding where your NCBI API Key is follow these directions: [https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us). Attempting to install this with conda did not end well so do it this way.

Use the script `download_venom_tsa_proteins.py` on the accession metadata list that already has split out for accessions that already have called proteins through the TSA with:
`scripts/download_venom_tsa_proteins.py metadata/SRA_TSA_venom_gland_accessions_protein_list.csv`. *For dumb reasons because of the way the Entrez tools have to be installed you must run it this way so that it calls python3 from /usr/bin/python3 and not python3 through an anaconda installation for example*.

The script will output all proteins in the TSA for each species in a fasta file named with the corresponding species name for each accession. For accessions that proteins couldn't be found, you can look at the resulting `download_status.log` to see which ones failed and use that for downloading the transcribed RNA sequences directly and performing downstream protein prediction. The `-d` parameter adds a time delay with a default of 180 seconds because if you ping NCBI too many times you will get your IP blocked, I usually use 60 seconds.
### Download transcribed sequences for accessions without proteins and get predicted proteins
This download script also takes advantage of the Entrez Direct toolkit from NCBI so install and follow directions as above. This takes a metadata file that has the columns accession, species_name, and contig so it can download all contigs for that species accession. The accession prefix won't work for efetch for FASTA files for some reason, so this script creates a list of all ranges of accesion + contig number up to the total number of contigs for that species and outputs to a FASTA for that species. The script is modified from the above downloading script for proteins but default is fasta and has modifications for the contig number addition. Use it with `scripts/download_venom_fasta_accessions.py metadata/SRA_TSA_venom_gland_accessions_fastarna_list.csv`.

For either downloading proteins or fastas I've found that NCBI has regularly had downtimes and for things that should download successfully there are sometimes HTTP/curl errors. Try again the next day.

### Run transdecoder on transcribed rna sequences for accessions without proteins
Use the bash script `scripts/run_transdecoder.sh` to process all fastas in a directory with transcribed rna sequences to obtain predicted ORFs/proteins. This worked on the majority of the fasta accessions, there were maybe 10-15 that failed for weird reasons (a couple were fasta formatting, and I'm not going to mess with it right now). This script runs `transdecoder` with default parameters, more parameters can be added to make ORF calling stricter, for now I just wanted to see how things worked.

### Rename protein locus tags and pool all proteins together
Put all the protein fasta files from accessions that originally had proteins and accessions where proteins were predicted from `transdecoder` into a single directory. Use the script `scripts/reorganize_fastas.py input_dir output.fasta` to reformat the headers to carry through the species name (which is the name of each fasta file) and pool all proteins together.

## Search against Pfam toxin domains and subsequently annotate with the entire Pfam database
This first analysis tested if we could quickly find candidate toxins among the pooled venom proteins by searching against Pfam toxin domains. Then we hypothesized that some toxin hits would have additional domains that are possibly the result of gene fusions and could act as facilitator proteins, so we asked this by annotating those hits with the entire Pfam database.

### Prepare toxin database
First we have some prep work for the Pfam database before we begin our searches.
1. Download the Pfam database from InterPro [https://www.ebi.ac.uk/interpro/download/pfam/](https://www.ebi.ac.uk/interpro/download/pfam/) and unzip
2. Create a conda environment with `hmmer` and dependencies: `conda create -n hmmer` activate it and install with `mamba install hmmer biopython` - you will want the biopython installation for later in the same environment
3. Obtain the list of specific Pfams you want, such as those corresponding to toxins in a TXT file
4. Since the Pfam specific number in `Pfam-A.hmm` is for example `PF19193.3` when you only have `PF19193` in the list, grep for the correct number with: `for line in $(cat toxin_list.txt); do grep -w $line Pfam-A.hmm; done > pfam_hits.txt | awk -F " " '{print $2}' > toxin_pfam_accessions.txt`. This file is in the `metadata/` directory
5. Obtain the specific toxin pfam accessions in the txt file list from the full `Pfam-A.hmm` database with `hmmfetch`. You can either set this up so each individual HMM is in a separate file or all together in one, which with about 480 of them it's fine if they are together for downstream running purposes. Do this with `hmmfetch -f Pfam-A.hmm ../metadata/toxin_pfam_accessions.txt > toxin_pfam_profiles.hmm`

### Perform HMM searches against known toxin HMMs and summarize the proteins that pass a certain threshold
Run the script `scripts/run_parse_toxin_hmms.py dbs/toxin_pfam_profiles.hmm results/all_toxin_proteins.fasta` to get a list of hits that pass a default evalue cutoff of 1e-50 against the toxin Pfams. The summary file lists for each hit protein the locus tag, length, cysteine content, and the Pfam accession/name it hit against.

### Get all additional Pfam annotations besides the toxin Pfams
Run the script `scripts/run_all_pfam.py --pfam_db dbs/Pfam-A.hmm results/toxin_protein_hits.fasta metadata/toxin_pfam_accessions.txt` to first run the entire Pfam HMM database against the FASTA of the toxin protein hits, and then it removes hits that match the the toxin Pfam accessions, since this will be redundant with the summary outfile in the previous step. Therefore the summary file contains any additional Pfam annotation for the hit proteins besides the toxin accessions.

### Summarize results
The notebook `toxin_facilitator_report.Rmd` was used to generate an HTML report that merged the summary files for the toxin hits and additional Pfam hits for proteins that had an additional hit. The notebook also plots protein length against cysteine content and colors points by the Pfam accession.

## Clustering curated Uniprot toxin proteins with venom proteins
Although the above Pfam `hmmsearch` analysis was successful and produced the expected data, the results are too general to move forward with decisive action. This is because the domains pulled for toxins from the Pfam HMMs can still hit pretty broad things without specificity, so it wasn't clear if the hits we pulled were actually toxins. This second approach uses `mmseqs` and full-length, curated toxin proteins from Uniprot to "search" against the venom proteins

### Prepare the Uniprot toxin protein database and venom proteins database
The curated set of Uniprot toxin proteins likely contain identical or similar homologs, and we want to cluster these to representative sequences so that we don't have multiple identical/similar clusters to compare against when querying against the venom proteins. First we will cluster this database of proteins into representative sequences with `mmseqs easy-cluster dbs/uniprot-toxin-proteins.fasta uniprot-toxin-clustering tmp`. This creates a file `uniprot-toxin-clustering_rep_seq.fasta` which are the representative sequences from clustering, which contains 1645 sequences, whereas the original uniprot toxin proteins totaled 7736.

Next make mmseqs databases of the representative uniprot toxin sequences and the target venom proteins:
```
mmseqs createdb all_venom_proteins.fasta venom_proteins_db
mmseqs createdb uniprot-toxin-clustering_rep_seq.fasta uniprot_toxin_reps_db
```

### Perform `mmseqs search`, parse results and produce summary statistics
Run `mmseqs search uniprot_toxin_reps_db venom_proteins_db resultDB tmp` where the first DB is the query DB and the second DB is the target DB. There will be several result files. To convert the results to tab-delimited in BLAST-formatted output style, use `mmseqs convertalis uniprot_toxin_reps_db venom_proteins_db resultDB uniprot-vs-venom-proteins-results.tsv`, which creates an output that looks like:

```
Q75WF2	Dendroaspis_jamesoni_kaimosae_GHPB01.transdecoder.pep_id_GHPB01000089.1.p1	0.429	353	200	0	5	357	39	389	1.182E-78	269
Q75WF2	Dendroaspis_jamesoni_jamesoni_GHPA01.transdecoder.pep_id_GHPA01000008.1.p1	0.445	332	184	0	26	357	18	349	1.441E-77	266
Q75WF2	Dendroaspis_polylepis_GHPD01.transdecoder.pep_id_GHPD01000149.1.p1	0.441	333	185	0	26	357	18	350	6.878E-77	264
Q75WF2	Bungarus_multicinctus_GIKH01.transdecoder.pep_id_GIKH01047790.1.p1	0.440	333	185	0	26	357	19	351	1.285E-76	263
Q75WF2	Naja_nigricollis_GIJG01.transdecoder.pep_id_GIJG01000137.1.p1	0.421	353	203	0	5	357	31	382	1.757E-76	263
Q75WF2	Naja_sumatrana_GIJN01.transdecoder.pep_id_GIJN01000117.1.p1	0.417	353	205	0	5	357	22	373	3.997E-75	259
Q75WF2	Ampulex_compressa_GFCP01.transdecoder.pep_id_GFCP01038013.1.p4	0.394	364	205	0	18	357	21	384	4.681E-71	247
Q75WF2	Ampulex_compressa_GFCP01.transdecoder.pep_id_GFCP01038014.1.p4	0.394	364	205	0	18	357	21	384	4.681E-71	247
Q75WF2	Ampulex_compressa_GFCP01.transdecoder.pep_id_GFCP01038017.1.p4	0.394	364	205	0	18	357	21	384	4.681E-71	247
Q75WF2	Ampulex_compressa_GFCP01.transdecoder.pep_id_GFCP01038018.1.p4	0.394	364	205	0	18	357	21	384	4.681E-71	247
```

From the mmseqs documentation:
> The file is formatted as a tab-separated list with 12 columns: (1,2) identifiers for query and target sequences/profiles, (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, (7-8, 9-10) domain start and end-position in query and in target, (11) E-value, and (12) bit score.

Now we have important information across multiple files - the resulting TSV from `mmseqs`, the original venom proteins fasta file that we can get protein length and cysteine content for the hits, and the uniprot metadata for each toxin accession. We can get a summary file containing for each hit the protein length, cysteine content, important hit information to each cluster (including alignment length, sequence identity, evalue, and bit score) along with the metadata for each accession (the long description) into one file with:

`python3 scripts/parse_mmseqs_results.py uniprot-vs-venom-proteins-results.tsv all_venom_proteins.fasta dbs/uniprot-toxin-proteins-metadata.tsv uniprot-vs-venom-proteins-summaries.tsv`.

This summary file is in the `results/2023-04-13-clustering-results` directory. The notebook `notebooks/toxin_clusters_summary_report.Rmd` was used to analyze the cluster results and show to the Strike team.

### Run against human proteins
The same procedure as above was run for a set of human proteins from Uniprot and searched against the representative toxin sequences by making an mmseqs database of the Uniprot human proteins and searching with default parameters. This will help to serve as a way to inform what thresholds we should be using for filtering bad hits and if there are things in the human proteome that are homologous to toxin proteins but don't quite serve a "toxic" function, such as hyaluronidases.

## Identification of protein length-outliers
The cluster results from the above mmseqs search are further refered to as
toxin-clusters'. Each cluster is characterized by a single toxin protein from the curated Uniprot database (refered to as the cluster reference toxin) and contains any venom protein that matched against the cluster reference toxin. Thus, we expect proteins from the same cluster to have the toxin domain and be of comparable size. Here, we want to identify length-outliers (but only longer proteins) in each cluster. The custom R script `Length-outliers-search.R` utilizes the file previously generated summary filr `uniprot-vs-all-venom-tick-proteins-summaries.tsv` and the curated Uniprot Venom proteins and Toxins database metadata file `Uniprot_Cluster_Info.csv` to generate the list of venom proteins that corresponds to our length-outlier definition (see below): `Venomproteins_ticks_toxins_outliers_20230428.csv`. The script relies on the R packages: `tiyverse`, `ggpubr`, and `plotly`.

To identify the proteins that are length-outliers:
1. Keep only one cluster hit per venom protein. If a venom protein was assigned to multiple clusters (ie had multiple hits in the curated Uniprot database), we only keep the hit that is associated with the lowest e-value.
2. Filter out any cluster that contains less than 5 venom proteins
3. Calculate or length evaluation metric: the protein length ratio. For each venom protein, its protein length ratio corresponds to the ratio of the length of the protein and the length of the cluster reference protein).
4. Identify any protein that meet the 'soft outliers' definition criteria: (i) length-ratio >=1.5 &, (ii) length-ratio > 1.5*IQ+Q3 ; where Q3 is the 3rd quartile of the distribution of the length-ratios of a cluster and IQ the interquartile range.
5. Visualize ratio distribution in each cluster and outliers. Plots are generated every 10 clusters (ordered according to the length of the cluster reference protein).

Protein outliers informations are then extracted and exported into the csv file: `Venomproteins_ticks_toxins_outliers_20230428.csv`.
Code for Figure3A (length-ratio distribution and outliers for 10 clusters) & Figure3B (heatmap of species and toxin were outliers are found) is part of the same custom R script. Figure 3B requires the extra file `Outliers_data_Figure1B.csv`. Figures can be found in the folder 'figures/'


## High-throughput extraction and annotation of accessory domains from outlier toxin protein hits
From exploring the clusters against the Uniprot toxin protein database, there are some protein hits that are much longer than the reference sequence or other protein hits in that cluster. These outliers were collected and listed in `results/2023-04-26-clustering-results/uniprot-toxin-proteins-outliers.xlsx`. This has the protein locus tag and what Uniprot toxin cluster it was assigned to. To extract the accessory domains:

1. Create a DIAMOND database of the Uniprot toxin proteins, and create a subset fasta file of the outlier proteins based on the locus tag from the concatenated input fasta file used for the initial clustering
2. Run `diamond blastp` of the outlier proteins against the Uniprot DIAMOND database
3. For each input protein, remove the part that hits against the Uniprot toxins, leaving just the accessory domains
4. Annotate these accessory domains against the Pfam HMMs, possibly also BLAST to NCBI-nr
5. Create clusters of the accessory domain proteins with `mmseqs` to get representative sequences and see if there is some conservation of accessory-domain types

### Create DIAMOND database of Uniprot toxin proteins
Create a DIAMOND conda environment with:

```
conda create -n diamond
conda activate diamond
mamba install diamond
```

And create a DIAMOND database of the Uniprot toxin proteins with: `diamond makedb --in dbs/uniprot-toxin-proteins.fasta --db uniprot-toxin-proteins.dmnd`.

Because this will include toxin proteins that are likely related homologs and the BLAST output will contain "duplicates" and be more difficult to parse, we can instead create a DIAMOND database of the representative sequences that were used for `mmseqs` clustering, which are in `results/2023-04-13-clustering/uniprot_clustering`. Make a DIAMOND database with `diamond makedb --in uniprot-toxin-clustering_rep_seq.fasta --db uniprot-representative-toxin-seqs.dmnd`
### Create subset fasta file of outlier proteins
The outlier proteins are listed in `results/2023-05-01-diamond-uniprot/lists/outlier-protein-hits-list.txt`, and I checked with `cat outlier-protein-hits-list.txt | sort | uniq` to make sure they are all unique locus tags, and there are a total of 702 protein locus tags in this list. Use the script `make-subset-fasta.py` to input a list of locus tags and create a new FASTA file with only those proteins, with `python3 scripts/make_subset_fasta.py lists/outlier-protein-hits-list.txt all_venom_tick_proteins.fasta subset_venom_tick_proteins_outliers.fasta`.

### Search subset proteins against the toxin DIAMOND database
Then with DIAMOND search the proteins against the toxin DIAMOND database with:
```
diamond blastp --db uniprot-toxin-proteins.dmnd \\
 --query subset_venom_tick_proteins_outliers.fasta \\
 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore corrected_bitscore qlen slen \\
 --out uniprot-toxins-vs-outlier-hit-proteins.txt
```

Perform a search against the representative Uniprot toxin sequences with:
```
diamond blastp --db uniprot-representative-toxin-seqs.dmnd \\
 --query subset_venom_tick_proteins_outliers.fasta \\
 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore corrected_bitscore qlen slen \\
 --out uniprot-toxins-rep-seqs-vs-outlier-hit-proteins.txt
 ```

### Remove toxin-like sequence from each hit protein to extract out "accessory-like" sequences
 There are a couple of possibilities that we could try to remove the toxin part of the protein to only have candidate accessory parts for each protein:

 1. Take the top hit based on e-value or bit score and extract the range of the query hit from that, and remove from the protein sequence
 2. Match up with the list Manon gave me to know what `mmseqs` cluster the protein ended up in, find the BLAST result for that uniprot protein hit, and extract the range of the query hit
 3. Take all possible blast hits and come up with the range from anything that encompasses those hits to remove everything and anything "toxin-like" to truly get hopefully accessory domains
 4. Another possibility that would create some redundancy but give more options is to create a multi-fasta file of each hit protein that has different "versions" of the remaining accessory domain based on what range each BLAST result for that protein looked like.

For an example of #4, this is what the first few BLAST hits for the same query protein look like against several Uniprot toxin DB sequences:

```
Dendroaspis_viridis_GHPC01.transdecoder.pep_id_GHPC01000014.1.p1	sp|Q98UF9|VM3H3_BOTJA	55.6	387	169	3	12	397	222	606	2.62e-155	446	521    	398	606	387
Dendroaspis_viridis_GHPC01.transdecoder.pep_id_GHPC01000014.1.p1	sp|P0C6E8|VM3G1_TRIGA	50.5	400	186	6	1	397	45	435	1.02e-134	387	453    	398	435	400
Dendroaspis_viridis_GHPC01.transdecoder.pep_id_GHPC01000014.1.p1	sp|P86092|VM3LB_BOTLC	50.6	326	152	4	66	389	1	319	1.68e-116	337	393    	398	324	326
Dendroaspis_viridis_GHPC01.transdecoder.pep_id_GHPC01000014.1.p1	sp|Q90WC0|VM2HS_GLOBR	48.9	272	134	3	1	271	44	311	9.84e-87	260	294
```

For the first one, the qlen is about half of the slen, where the slen is much longer. For the third hit, the query and db proteins are about the same length. So for these you would basically remove the entire query protein, there would be some leftover here and there. But for hits with a lower evalue for example:

```
Dendroaspis_viridis_GHPC01.transdecoder.pep_id_GHPC01000014.1.p1	sp|P31989|VM2JC_BOTJA	50.0	158	74	3	115	271	2	155	3.30e-48	155	167    	398	161	158
Dendroaspis_viridis_GHPC01.transdecoder.pep_id_GHPC01000014.1.p1	sp|Q0NZX7|VM3B4_BOTJA	67.0	97	32	0	189	285	1	97	1.24e-44	144	156    	398	97	97
Dendroaspis_viridis_GHPC01.transdecoder.pep_id_GHPC01000014.1.p1	sp|P0DJ43|VM3_MICIK	36.5	252	41	4	146	396	34	167	6.06e-43	142	149    	398	169	252
```

Where the query is about 2x or more longer than the db protein, there is quite a bit leftover that could be "accessory". These accessory parts could already be in the other DB proteins that are quite long, we would just have to do another BLAST against the DB proteins with the accessory parts subset out. This could potentially create redundant or duplicate "accessory" proteins in the FASTA file, but we could also cluster them into representative sequences using `mmseqs` like we originally did for the Uniprot toxin DB sequences in the first place (filtering for things above a certain length).

To test this I subset the BLAST output file down to just the locus tag, start, and stop columns with `awk -F "\t" '{print $1"\t"$7"\t"$8}' uniprot-toxins-vs-outlier-hit-proteins.txt > uniprot-toxins-outliers-ranges.txt`. Then I used the script `scripts/remove-blast-hit-from-protein.py` to take each line from the BLAST results, find the protein, remove that range of the sequence, and write the remaining sequence to a new FASTA file. I called the script with `python3 scripts/remove-blast-hit-from-protein.py uniprot-toxins-outliers-ranges.txt subset_venom_tick_proteins_outliers.fasta accessory_toxin_sequences.fasta`.

I spot checked about 10 random proteins to make sure it removed the correct part as expected. I only made a check to remove records that are empty and contain no amino acids (BLAST results where it hit the entire range of the protein) and kept regardless of length. Therefore this keeps really short proteins/peptides that are a few AAs, and long proteins that could be comprised of what was on either end of the toxin hit, not just before or after. I also added a check to remove identical proteins so there isn't a bunch of duplicates, although the next step would take care of that as well.

### Create representative sequences of these "accessory" domains
The previous step has a check ot remove identical protein sequences from being written to the multi-FASTA file of "accessory" parts of each protein that are left over from the toxin-like hit from the BLAST searches to the Uniprot toxin proteins. Now I'm going to cluster these into representative sequences since although identical sequences were removed, just from spot checking there are some that are still highly similar, and we want to simplify things here.

We will create representativ sequences with the `mmseqs easy-cluster` program like we did for creating representative sequences of the Uniprot toxin proteins with `mmseqs easy-cluster accessory_nondup_toxin_sequences.fasta accessory-toxin-proteins-clustering tmp`. The original fasta file of "accessory" file had 2068 sequences, and through clustering to representative sequences this new FASTA file has 377 sequences in the file `accessory-toxin-proteins-clustering_rep_seq.fasta` and you can see what ended up in what cluster in `accessory-toxin-proteins-clustering_cluster.tsv`.

### Annotate the representative "accessory" proteins with the Pfam database
Now I will take the representative 377 proteins that are "accessory" from the toxins and annotate against the Pfam.hmm database. We can use the `scripts/run_parse_hmms.py` script that we used for the original part of this project with:

```
python3 run_parse_hmms.py \
    dbs/Pfam-A.hmm \
    accessory-toxin-proteins-clustering_rep_seq.fasta \
    --evalue 1e0
    --output_file accessory-rep-seqs-v-pfam.fasta
```

## Accessory sequences analysis
We further filter/investigate the accessory sequences clusters identified in `accessory-toxin-proteins-clustering_cluster.tsv`. The objective is to characterize the diversity of species and toxin within each clusters. In other word we want to understand if accessory sequences in each cluster are associated with the same species and if they are accessory sequences of the same toxin (as identified by the toxin-clusering).
This analysis is performed in the custom R script `Accessory_sequence_diversity_analysis.R` and requires the input files: ` Clusters_Accessory_05022023.csv` (reformated output of the accessory sequences clustering), `Cluster_Acc_Ref_PFAM.csv` (Pfam annotations of the accessory sequences clusters) and `Origin_data_accessory_seq_clusters.csv ` (Toxin information for accessory sequences - indicates for each present accessory sequence what is the toxin an accessory sequence is expected to be associated with according to the cluster reference toxin of the toxin-cluster where the accessory sequence is from). This script produces the Pub's Figures 4 and 5 as well as Table 2.
This scripts uses the R packages: `tidyverse`, `ggplot2`,`dplyr` and `plotly`.

### Cluster filtering based on Pfam annotation
Toxins span a broad range of protein family. Pfam consistently associated with toxins (toxin-associated Pfam) can be inferred from the curated Uniprot Venom protein and toxins database and compiled with the Pfam list identified in the following work: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5655178/ . According to this list: `list`, we identify any accessory sequence clusters whose representative sequence has been annotated with a toxin-associated Pfam and filter them out.
Code for Figure 4 (histogram of the accessory sequences clusters size colored based on the Pfam annotation status) is present in this section


### Species and toxin diversity analysis
In each remaining cluster, we further identify: from how many different outlier sequences the accessory sequences have been extracted and visualize the distribution of the number of outlier/cluster in Figure 5A.
For the rest of the diversity analysis, we filter out any accessory sequences that is assocaited with a single outlier. Then, we continue the analysis by calculating and plotting
    - the number of species represented / cluster (Fig 5B)
    - the number of different toxins (based on toxin-clusters reference toxin) associated with the accessory sequences / cluster (Fig 5C)
