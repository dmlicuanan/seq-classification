#!/bin/bash

# this script is for downloading GenBank flat files or fasta files from NCBI using Entrez Direct (https://www.ncbi.nlm.nih.gov/books/NBK179288/)
# gene: COI sequences from nuccore database
# species: marine eukaryotes 

# species list is in the file .../OneDrive/2022_kraken/input/wormstaxlist
# wormstaxlist was derived from 2022-07-01 copy of World Register of Marine Species
# it includes species under kingdoms Fungi, Animalia, Chromista, Protozoa, and Plantae
# it is restricted to species with "accepted" taxonomicStatus





########## version 1 -- last edit 22/07/2022

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





########## draft for version 2 

# set working directory (OneDrive)
# cd /mnt/c/Users/Ardea\ Licuanan/OneDrive/2022_kraken/
# working directory for trial 
cd /mnt/d/Documents/NGS/entrez/tmp/ncbi_trial

# use API key of amlicuanan@alum.up.edu.ph = 81af93ade140a45a355b621d22a8692a5408
export NCBI_API_KEY=81af93ade140a45a355b621d22a8692a5408

# create directory where GenBank files will be saved:
mkdir wormstaxlist_gb

# copy wormstaxlist to temporary folder for trial 
cp /mnt/d/Documents/NGS/entrez/wormstaxlist wormstaxlist

# check taxa with special characters 
# all special characters must be URL encoded -- https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
grep -v "^[a-zA-Z -.]*$" wormstaxlist 

# remove strings [sensu lato] and (incertae sedis) + characters × and ,
# delete lines with [non-Uristidae]
# convert all characters to ascii
# remove single and double quotes
sed -e "s/\( \[sensu lato\]\| (incertae sedis)\|[×,]\)//" wormstaxlist | sed '/\[non-Uristidae\]/d' | iconv -f utf-8 -t ascii//translit | tr -d \'\" > wormstaxlist_cleaned
# check if all special characters are removed
grep -v "^[a-zA-Z -.]*$" wormstaxlist_cleaned 

# add FIELDS to each line in taxa list
sed 's/$/ [ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])/' wormstaxlist_cleaned > wormstaxlist_fin
head -n 500 wormstaxlist_fin > subset

# alternative 1 to original:
# time parallel -j 4 --delay 0.34 --keep-order 'esearch -db nuccore -query {} | xtract -pattern Count -element Count | xargs -I [] echho "[]_{}"  >> wormstaxlist_counts' :::: ./input/wormstaxlist_fin

# prints number of hits:
# echo "a" | xargs -I [] echo "[]_b"
# parallel -j 4 --delay 0.34 --keep-order 'esearch -db nuccore -query {} | \
# xtract -pattern Count -element Count | \
# xargs -I [] echo "[] {}" >> counts' :::: subset

# alternative 2 to original:
# function outputs hits and query in file
# counts() {
	# query="$1"
	# IFS=$'\n'
	# for i in $(esearch -db nuccore -query "$query" | xtract -pattern Count -element Count)
		# do echo "${i}, $query" >> counts
		# done
	# }
# export -f counts

# run function in parallel (10 jobs at a time, delayed start of 1.1 seconds)
# time parallel -j 10 --delay 1.1 counts :::: subset

# function for downloading GB flat file
mydownload() {
	query="$1"
	IFS=$'\n'
	search=$(esearch -db nuccore -query "$query") 
	count=$(echo "$search" | xtract -pattern Count -element Count)
	if [ $count -gt 0 ]
	then 
	filename=$(echo "$query" | sed 's/ /_/;s/\[.*//;s/ //')
	echo "$search" | efetch -format gb > ./wormstaxlist_gb/${filename}_ncbi.gb
	fi	
	}
export -f mydownload

# run in parallel (10 jobs at a time, delayed start of 1.1 seconds)
# save errors to log_ncbi.txt
time parallel -j 10 --delay 1.1 mydownload :::: subset 2>&1 | tee -a log_ncbi.txt




########## version 2 -- last edit 29/07/2022

# set working directory (OneDrive)
cd /mnt/c/Users/Ardea\ Licuanan/OneDrive/2022_kraken/

# use API key of amlicuanan@alum.up.edu.ph 
export NCBI_API_KEY=81af93ade140a45a355b621d22a8692a5408
# other API keys:
# 7d21c9bc180e2c7f4c118b010e6f236f3008 
# f4aef0d264bef47873c8c44329570bda3d08 
# 532b235a5e1021b4e209f1d22f6b69258e08 
# bfeb59eb00909bc3522de6eb2214bd769508 
# beccc4a5d3d81297a8625322e2cce745a508 
# c0f0673298193bfbb7dfacb8b362375a2008 

# create directory where GenBank files will be saved:
mkdir wormstaxlist_gb

# check taxa with special characters 
# all special characters must be URL encoded -- https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
grep -v "^[a-zA-Z -.]*$" ./input/wormstaxlist 

# remove strings [sensu lato] and (incertae sedis) + characters × and ,
# delete lines with [non-Uristidae]
# convert all characters to ascii
# remove single and double quotes
sed -e "s/\( \[sensu lato\]\| (incertae sedis)\|[×,]\)//" ./input/wormstaxlist | sed '/\[non-Uristidae\]/d' | iconv -f utf-8 -t ascii//translit | tr -d \'\" > wormstaxlist_cleaned
# check if all special characters are removed
grep -v "^[a-zA-Z -.]*$" wormstaxlist_cleaned 

# add FIELDS to each line in taxa list
sed 's/$/ [ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])/' wormstaxlist_cleaned > wormstaxlist_fin

# create function for downloading GB flat file
mydownload() {
	query="$1"
	IFS=$'\n'
	search=$(esearch -db nuccore -query "$query") 
	count=$(echo "$search" | xtract -pattern Count -element Count)
	if [ $count -gt 0 ]
	then 
	filename=$(echo "$query" | sed 's/ /_/;s/\[.*//;s/ //')
	echo "$search" | efetch -format gb > ./wormstaxlist_gb/${filename}_ncbi.gb
	fi	
	}
export -f mydownload

# run in parallel (10 jobs at a time, delayed start of 1.1 seconds)
# save errors to log_ncbi.txt
time parallel -j 10 --delay 1.1 mydownload :::: wormstaxlist_fin 2>&1 | tee -a log_ncbi.txt