#!/bin/bash

# this script is for downloading GenBank flat files or fasta files from NCBI using Entrez Direct (https://www.ncbi.nlm.nih.gov/books/NBK179288/)
# gene: COI sequences from nuccore database
# species: marine eukaryotes 

# species list is in the file .../OneDrive/2022_kraken/input/wormstaxlist
# wormstaxlist was derived from 2022-07-01 copy of World Register of Marine Species
# it includes species under kingdoms Fungi, Animalia, Chromista, Protozoa, and Plantae
# it is restricted to species with "accepted" taxonomicStatus

# set working directory (OneDrive)
cd /mnt/c/Users/Ardea\ Licuanan/OneDrive/2022_kraken/

# create directory where GenBank files will be saved:
mkdir wormstaxlist_gb
# create directory where fasta files will be saved:
mkdir wormstaxlist_fasta

# add FIELDS to each line in wormstaxlist
sed 's/$/ [ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])/' ./input/wormstaxlist > ./input/wormstaxlist_fin

# record which species have hits 
time parallel -j 4 --delay 0.34 --keep-order 'esearch -db nuccore -query {} | xtract -pattern Count -element Count >> wormstaxlist_counts' :::: ./input/wormstaxlist_fin
# note: for 1000 species, this took real-time: 84m49.175s
# --delay is for preventing warnings about too many requests

# print only species which have hits; we use this list later for efetch
paste wormstaxlist_counts ./input/wormstaxlist_fin | awk '$1 > 0' | cut -f 2- > ./input/wormstaxlist_whits

# download GenBank flat files for each species
time parallel -j 4 --delay 0.34 'esearch -db nuccore -query {} | efetch -format gb > ./wormstaxlist_gb/{=s/ /_/;s/\[.*//;s/ //;=}.gb' :::: ./input/wormstaxlist_whits

# download fasta file for each species
time parallel -j 4 --delay 0.34 'esearch -db nuccore -query {} | efetch -format fasta > ./wormstaxlist_fasta/{=s/ /_/;s/\[.*//;s/ //;=}.fas' :::: ./input/wormstaxlist_whits
# note: for 455 species, this took realtime: 20m15.668s
