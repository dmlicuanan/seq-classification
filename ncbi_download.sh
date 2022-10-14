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




########## version 2 for running in Sir's machine -- last edit 29/07/2022

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




########## version for "back-up run" from 2022/08/01 to 2022/08/06
# set working directory 
cd /mnt/d/Documents/NGS/entrez/

# use API key of amlicuanan@alum.up.edu.ph 
export NCBI_API_KEY=81af93ade140a45a355b621d22a8692a5408

# create directory where GenBank files will be saved:
# mkdir wormstaxlist_gb

# check taxa with special characters 
# all special characters must be URL encoded -- https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
# grep -v "^[a-zA-Z -.]*$" wormstaxlist 

# remove strings [sensu lato] and (incertae sedis) + characters × and ,
# delete lines with [non-Uristidae]
# convert all characters to ascii
# remove single and double quotes
# sed -e "s/\( \[sensu lato\]\| (incertae sedis)\|[×,]\)//" wormstaxlist | sed '/\[non-Uristidae\]/d' | iconv -f utf-8 -t ascii//translit | tr -d \'\" > wormstaxlist_cleaned
# check if all special characters are removed
# grep -v "^[a-zA-Z -.]*$" wormstaxlist_cleaned 

# add FIELDS to each line in taxa list
# sed 's/$/ [ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])/' wormstaxlist_cleaned > wormstaxlist_fin

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

# create subsets of 10,000 and put in separate directory
# mkdir wormstaxlist_gb/subsets
# split --verbose -l10000 wormstaxlist_fin ./wormstaxlist_gb/subsets/subset.

# run in parallel (10 jobs at a time, delayed start of 1.1 seconds)
# save errors to log_ncbi.txt appended by subset identifier
# change subset files manually
time parallel -j 10 --delay 1.1 mydownload :::: ./wormstaxlist_gb/subsets/subset.ab 2>&1 | tee -a ./wormstaxlist_gb/subsets/log_ncbi_ab.txt 
# 8/1/2022
# subset.aa 10:04 AM (186m34.434s)
# subset.ab 1:27 PM (95m51.774s) only 1612 downloaded; re-run on 8/8/2022 6:04 PM (197m39.768s) 2782 downloaded after re-runs
# subset.aw 6:26 PM (45m14.437s)
# subset.ac 7:19 PM (195m9.809s)
# subset.ad 10:42 PM (195m8.691s)
# 8/2/2022
# subset.ae 9:19 AM (207m3.172s)
# subset.af 5:28 PM (194m25.961s)
# subset.ag 8:54 PM (195m13.427s)
# 8/3/2022
# subset.ah 9:22 AM (195m59.597s)
# subset.ai 5:16 PM (194m30.055s)
# subset.aj 8:53 PM (194m25.263s)
# subset.ak 12:14 AM (193m20.197s)
# 8/4/2022
# subset.al 9:48 AM (194m31.344s)
# subset.am 1:11 PM (194m51.515s)
# subset.an 5:49 PM (193m40.201s)
# subset.ao 9:40 PM (193m18.177s)
# 8/5/2022
# subset.ap 9:12 AM (194m50.979s)
# subset.aq 12:40 PM (194m52.813s)
# subset.ar 5:40 PM (193m2.918s)
# subset.as 10:10 (192m59.591s)
# 8/6/2022
# subset.at 8:13 AM (193m27.705s)
# subset.au 12:15 PM (193m13.827s)
# subset.av 3:57 PM (193m15.131s)

# 2022/08/08: check log files for species that need to be re-run
cd /mnt/d/Documents/NGS/entrez/wormstaxlist_gb/subsets/
# concatenate log files into a single file
cat *.txt >> log_ncbi_compiled.txt

# match string between -term and [ORGN] in log files
# replace underscore with space and append other fields
# https://stackoverflow.com/questions/13242469/how-to-use-sed-grep-to-extract-text-between-two-words
grep -o -P '(?<= -term ").*?(?= \[ORGN\])' log_ncbi_compiled.txt | sed 's/_/ /g' | sed 's/$/ [ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])/' | sort -u > /mnt/d/Documents/NGS/entrez/wormstaxlist_rerun
# re-run for species in log_ncbi_ab.txt (after 2022/08/08 re-run)
grep -o -P '(?<= -term ").*?(?= \[ORGN\])' log_ncbi_ab.txt | sed 's/_/ /g' | sed 's/$/ [ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])/' | sort -u > /mnt/d/Documents/NGS/entrez/wormstaxlist_rerun_ab

# perform re-run -- 8/8/2022 (14m31.417s); 
cd /mnt/d/Documents/NGS/entrez
parallel -j 10 --delay 1.1 mydownload :::: wormstaxlist_rerun 2>&1 | tee -a ./wormstaxlist_gb/subsets/log_ncbi_rerun.txt 
# re-run for species in log_ncbi_ab.txt (after 2022/08/08 re-run)
# parallel -j 10 --delay 1.1 mydownload :::: wormstaxlist_rerun_ab 2>&1 | tee -a ./wormstaxlist_gb/subsets/log_ncbi_rerun.txt 
# no errors in log file, so will delete that and wormstaxlist_rerun_ab
# rm log_ncbi_rerun.txt /mnt/d/Documents/NGS/entrez/wormstaxlist_rerun_ab

# 149 of 150 files are existing and need to be overwritten in their respective folders
cd /mnt/d/Documents/NGS/entrez/wormstaxlist_gb
# move all re-run species from log files to their respective folders
for i in $(ls *_ncbi.gb) 
do
spec=$(echo "$i" | sed 's/_/ /;s/_ncbi.gb//') 
path=$(grep -w "$spec" /mnt/d/Documents/NGS/entrez/wormstaxlist_gb/subsets/subset* | cut -d : -f 1 | sed 's/subsets\///')
mv "${i}" "${path}"
done 

# check if last files downloaded per subset reached ~10K index
cd /mnt/d/Documents/NGS/entrez/wormstaxlist_gb/subsets
grep "Terschellingia longicaudata" -n subset.aa
grep "Metschnikowia tuberculata" -n subset.ab
grep "Diplasterias brucei" -n subset.ac
grep "Pervagor melanocephalus" -n subset.ad
grep "Diuronotus aspetos" -n subset.ae
grep "Abudefduf hoefleri" -n subset.af
grep "Euphysa japonica" -n subset.ag
grep "Speleobregma lanzaroteum" -n subset.ah
grep "Ostrincola japonica" -n subset.ai
grep "Macroparalepis johnfitchi" -n subset.aj
grep "Ancistrum crassum" -n subset.ak
grep "Euphilomedes morini" -n subset.al
grep "Glabella rosadoi" -n subset.am
grep "Timoclea scabra" -n subset.an
grep "Marionia cucullata" -n subset.ao
grep "Homoieurete macquariense" -n subset.ap
grep "Cetoscarus ocellatus" -n subset.aq
grep "Sermyla riquetii" -n subset.ar
grep "Amathia tertia" -n subset.as
grep "Philypnodon grandiceps" -n subset.at
grep "Ophryotrocha urbis" -n subset.au
grep "Spongilla manconiae" -n subset.av
grep "Hyalinoecia longibranchiata" -n subset.aw
# except for subset.ab whose last download was for the 5116th species Metschnikowia tuberculata, all reached expected index

# check that there are no hits for species beyond index 5116 in subset.ab; delete subset.ab_check later 
tail -n 4884 subset.ab > subset.ab_check
cd /mnt/d/Documents/NGS/entrez/wormstaxlist_gb/subsets
parallel -j 10 --delay 1.1 --keep-order 'esearch -db nuccore -query {} | xtract -pattern Count -element Count | xargs -I [] echo "[] {}" >> counts_ab_check' :::: subset.ab_check
# check how many species with hits were missed in first run (1162 species)
awk '$1 > 0' counts_ab_check | wc -l

# need confirmation that the rest of species (222,316-31,689=190,627) in wormstaxlist_fin have 0 hits
cd /mnt/d/Documents/NGS/entrez/wormstaxlist_gb/
# create list of species with hits, then exclude those from wormstaxlist_fin
find ./subset.*/ -name "*_ncbi.gb" | cut -d \/ -f 3 | sed 's/_/ /;s/_ncbi.gb//;s/$/ [ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])/' > /mnt/d/Documents/NGS/entrez/wormstaxlist_fin_whits

# create list of species with hits, then exclude those from wormstaxlist_fin
find ./subset.*/ -name "*_ncbi.gb" | cut -d \/ -f 3 | sed 's/_/ /;s/_ncbi.gb//' > /mnt/d/Documents/NGS/entrez/wormstaxlist_fin_whits

# create list of species that supposedly have no hits
cd /mnt/d/Documents/NGS/entrez/
# using whole words (-w) from file (-f) wormstaxlist_fin_whits, select non-matching lines (-v)
grep -v -w -f wormstaxlist_fin_whits wormstaxlist_fin > wormstaxlist_fin_for_counts

# create function for checking hits per species
hitcount() {
query="$1"
IFS=$'\n'
count=$(esearch -db nuccore -query "$query" | xtract -pattern Count -element Count) 
spec=$(echo "$query" | cut -d \[ -f 1)
echo "$count $spec" >> counts_wormstaxlist_fin_for_counts
}
export -f hitcount

# initial run started 2:08 PM 2022/08/09
export NCBI_API_KEY=81af93ade140a45a355b621d22a8692a5408
parallel -j 10 --delay 1 --keep-order hitcount :::: wormstaxlist_fin_for_counts 2>&1 | tee -a log_counts.txt

# there are duplicates in wormstaxlist_fin_for_counts (190,627); we remove those without sorting:
cat -n wormstaxlist_fin_for_counts | sort -u -k 2 | sort -n | cut -f 2- > wormstaxlist_fin_for_counts_uniq
# there are 190,388 unique entries

# run codes below for all successive runs:
# reduce unique species in wormstaxlist_fin_for_counts to remove those which have already have counts
# remove count in counts_wormstaxlist_fin_for_counts
# append species with [ORGN] to be more specific in grep
# get list of (unique) species to input as file of fixed-strings in (inverted) grep
# <() is called process substitution -- provides a way to pass the output of a command to another command when using a pipe is not possible
grep -v -w -F -f <(cut -d ' ' -f 2- counts_wormstaxlist_fin_for_counts | sed 's/$/\[ORGN\]/' | sort -u) wormstaxlist_fin_for_counts_uniq > wormstaxlist_fin_for_counts_reduced
# execute hitcount function:
parallel -j 10 --delay 1 --keep-order hitcount :::: wormstaxlist_fin_for_counts_reduced 2>&1 | tee -a log_counts.txt 
# number of lines of counts_wormstaxlist_fin_for_counts after initial run = 190389

# after initial run for counts, print species with hits > 1
awk '$1 > 0' counts_wormstaxlist_fin_for_counts 
# there are species which do not have counts yet but were printed in the file
# despite having no counts, they are still excluded from wormstaxlist_fin_for_counts_reduced 
awk '$1 > 0' counts_wormstaxlist_fin_for_counts | cut -d ' ' -f 2-
# these species need to be extracted from counts_wormstaxlist_fin_for_counts for rerun
# select lines that start with space and remove leading spaces and trailing spaces
awk '$1 > 0' counts_wormstaxlist_fin_for_counts | grep "^ " | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' 
# merge with list of species that occurred in log_counts.txt
# use sort -u to get unique values then append with fields
cat <(awk '$1 > 0' counts_wormstaxlist_fin_for_counts | grep "^ " | sed 's/^[[:space:]]*//;s/[[:space:]]*$//') <(grep -o -P '(?<= -term ").*?(?= \[ORGN\])' log_counts.txt | sed 's/_/ /g') | sort -u | sed 's/$/ [ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])/' > wormstaxlist_fin_for_counts_rerun
# lines of wormstaxlist_fin_for_counts_rerun = 1068

# run hit counter again for species with blank counts and species in log_counts.txt
# hitcount function will append new counts to counts_wormstaxlist_fin_for_counts so there will be duplicate species (some with blank counts, some with counts)
parallel -j 10 --delay 1 --keep-order hitcount :::: wormstaxlist_fin_for_counts_rerun 2>&1 | tee -a log_counts_rerun.txt 
# note: if there are errors in log_counts_rerun.txt, need to rerun those
# delete log_counts_rerun.txt if empty

# check that species originally with no count, now have counts
# use comm to compare two sorted files and only print species unique to file 1
# file 1 is list of species with empty counts
# file 2 is list of species with counts
comm <(grep -v '^[0-9]' counts_wormstaxlist_fin_for_counts | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | sort -u) <(grep '^[0-9]' counts_wormstaxlist_fin_for_counts | cut -d ' ' -f 2- | sed 's/[[:space:]]*$//' | sort -u) -23
# result should be empty

# check that all species that need their counts checked have counts
# file 1 is the list of species run for counts minus the fields
# file 2 is list of species with counts
comm <(sed 's/ \[ORGN.*//' wormstaxlist_fin_for_counts | sort -u) <(grep '^[0-9]' counts_wormstaxlist_fin_for_counts | cut -d ' ' -f 2- | sed 's/[[:space:]]*$//' | sort -u) -23
# result should also be empty

# 2022/08/11
# download GenBank files of species in counts_wormstaxlist_fin_for_counts which have counts
# these were missed in the initial (subset) download
# select species in list which have counts > 0 and append with fields
grep '^[1-9]' counts_wormstaxlist_fin_for_counts | cut -d ' ' -f 2- | sed 's/$/[ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])/' | sort -u > wormstaxlist_fin_whits_missed
# use mydownload function again to download GBs
cd /mnt/d/Documents/NGS/entrez/
export NCBI_API_KEY=81af93ade140a45a355b621d22a8692a5408
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
# run download
# delete log_ncbi_missed.txt if empty after run
parallel -j 10 --delay 1.1 mydownload :::: wormstaxlist_fin_whits_missed 2>&1 | tee -a log_ncbi_missed.txt 

# for-loop for moving .gb files from wormstaxlist_gb to their respective folders operates on the assumption that a species belongs to only one subset folder
# however, there are duplicates (and triplicates) in wormstaxlist_fin, meaning one species may be listed > 1 subset e.g.:
for i in $(sort wormstaxlist_fin_whits | uniq -d | sed 's/ /_/;s/$/_ncbi.gb/'); do find ./wormstaxlist_gb/subset.*/ -name "${i}"; done
# need to check that .gb files to be moved are not part of duplicate/triplicate list
# use comm
# file 1: species duplicated/triplicated in wormstaxlist_fin
# file 2: species with .gbs that need to be moved
comm <(sed 's/ \[ORGN.*//' wormstaxlist_fin | sort | uniq -d) <(find ./wormstaxlist_gb -maxdepth 1 -name "*_ncbi.gb" | cut -d \/ -f 3 | sed 's/_/ /;s/_ncbi.gb//' | sort -u) -12
# common species are printed; proceed only if empty result

# number of species before moving = 31689
find ./wormstaxlist_gb/subset.*/ -name "*_ncbi.gb" | wc -l
# number of species to be moved 
find ./wormstaxlist_gb -maxdepth 1 -name "*_ncbi.gb" | wc -l

# move .gbs of missed species to their respective folders
cd /mnt/d/Documents/NGS/entrez/wormstaxlist_gb
# copy of earlier code for moving files
for i in $(ls *_ncbi.gb) 
do
spec=$(echo "$i" | sed 's/_/ /;s/_ncbi.gb//') 
path=$(grep -w "$spec" /mnt/d/Documents/NGS/entrez/wormstaxlist_gb/subsets/subset* | cut -d : -f 1 | sed 's/subsets\///')
mv "${i}" "${path}"
done 

# check how many .gb files there are = 31721
find ./subset.*/ -name "*_ncbi.gb" | cut -d \/ -f 3 | wc -l 
# check how many (unique) .gb files there are = 31708
find ./subset.*/ -name "*_ncbi.gb" | cut -d \/ -f 3 | sort -u | wc -l

# check if duplicates are identical 
# https://unix.stackexchange.com/questions/236581/how-to-pipe-output-of-find-as-input-for-diff
cd /mnt/d/Documents/NGS/entrez/
for i in $(comm <(find ./wormstaxlist_gb/subset.*/ -name "*_ncbi.gb" | cut -d \/ -f 4 | sort -u) <(sort wormstaxlist_fin | uniq -d | sed 's/ \[ORGN.*/_ncbi.gb/;s/ /_/') -12)
do
diff -q $(find ./wormstaxlist_gb/subset.*/ -name "${i}")
done
# will print report if files differ; -s another alternative

# create directory where all _ncbi.gb will be copied to for compilation
# since duplicate files are identical, overwriting will be okay
cd /mnt/d/Documents/NGS/entrez/wormstaxlist_gb
mkdir compiled_gb
cp $(find ./subset.*/ -name "*_ncbi.gb") compiled_gb



# 2022/08/15
# re-download badfiles discovered in gb_parse (listed in wormstaxlist_badfiles)
# use mydownload function again to download GBs
cd /mnt/d/Documents/NGS/entrez/
export NCBI_API_KEY=81af93ade140a45a355b621d22a8692a5408
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
# run download (modified)
parallel -j 2 --delay 2 mydownload :::: wormstaxlist_badfiles

# copy new .gbs of "bad files" to their respective folders
# each file is not part of duplicated species in wormstaxlist_fin
cd /mnt/d/Documents/NGS/entrez/wormstaxlist_gb
# modification of earlier code for moving files 
# we will only copy them to their respective subset folders then move them to compiled_gb
for i in $(ls *_ncbi.gb) 
do
spec=$(echo "$i" | sed 's/_/ /;s/_ncbi.gb//') 
path=$(grep -w "$spec" /mnt/d/Documents/NGS/entrez/wormstaxlist_gb/subsets/subset* | cut -d : -f 1 | sed 's/subsets\///')
cp "${i}" "${path}"
done 
# move files to compiled_gb
mv *_ncbi.gb compiled_gb




########## download of 16S sequences (2022/10/14)
# set working directory 
cd /mnt/d/Documents/NGS/entrez/

# use API key of amlicuanan@alum.up.edu.ph 
export NCBI_API_KEY=81af93ade140a45a355b621d22a8692a5408

# create directory where GenBank files will be saved:
mkdir -p wormstaxlist_gb_16S/raw

# check taxa with special characters 
# all special characters must be URL encoded -- https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
# grep -v "^[a-zA-Z -.]*$" wormstaxlist 

# remove strings [sensu lato] and (incertae sedis) + characters × and ,
# delete lines with [non-Uristidae]
# convert all characters to ascii
# remove single and double quotes
# sed -e "s/\( \[sensu lato\]\| (incertae sedis)\|[×,]\)//" wormstaxlist | sed '/\[non-Uristidae\]/d' | iconv -f utf-8 -t ascii//translit | tr -d \'\" > wormstaxlist_cleaned
# check if all special characters are removed
# grep -v "^[a-zA-Z -.]*$" wormstaxlist_cleaned 

# add FIELDS to each line in taxa list
sed 's/$/ [ORGN] AND (16S [GENE] OR 16s [GENE] OR 16S ribosomal RNA [GENE] OR 16S rRNA [GENE])/' wormstaxlist_cleaned > wormstaxlist_fin_16S

# set variable
export gbout=/mnt/d/Documents/NGS/entrez/wormstaxlist_gb_16S

# create function for downloading GB flat file
mydownload() {
query="$1"
IFS=$'\n'
search=$(esearch -db nuccore -query "$query") 
count=$(echo "$search" | xtract -pattern Count -element Count)
if [ $count -gt 0 ]
then 
filename=$(echo "$query" | sed 's/ /_/;s/\[.*//;s/ //')
echo "$search" | efetch -format gb > $gbout/raw/${filename}_ncbi.gb
fi	
echo "$query" >> $gbout/wormstaxlist_fin_16S_done
}
export -f mydownload

# initial run:
# run in parallel (10 jobs at a time, delayed start of 1.1 seconds)
# save errors to log_ncbi_16S.txt 
parallel -j 10 --delay 1.1 mydownload :::: wormstaxlist_fin_16S 2>&1 | tee -a $gbout/log_ncbi_16S.txt 

# for succeeding runs, remove species that have already been run
comm <(sort -u wormstaxlist_fin_16S) <(sort -u $gbout/wormstaxlist_fin_16S_done) > wormstaxlist_fin_16S_reduced -23
parallel -j 10 --delay 1.1 mydownload :::: wormstaxlist_fin_16S_reduced 2>&1 | tee -a $gbout/log_ncbi_16S.txt 