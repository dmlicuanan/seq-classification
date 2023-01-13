# export conda environments for future use:
cd /mnt/d/Documents/NGS/conda_envs
# ngs
conda activate ngs
conda env export > ngs.yml
# py3.7
conda activate py3.7
conda env export > py3.7.yml
# to recreate environment:
conda env create -f ngs.yml
conda env create -f py3.7.yml

# create directories for outputs
mkdir /mnt/d/Documents/NGS/out/Manila_Bay
# create subdirectories for outputs
cd /mnt/d/Documents/NGS/out/Manila_Bay
mkdir fastqc multiqc transients bbmerge bbduk fastqscreen kraken

# set shell variables
# out = directory of outputs
out=/mnt/d/Documents/NGS/out/Manila_Bay
# in = directory of raw reads
in=/mnt/d/Documents/NGS/input/fastqFiles/Manila_Bay
# directory of bbmapresources
bbmapResource=/home/dmlicuanan/miniconda3/envs/ngs/opt/bbmap-38.18/resources
# directory of transient files
transients=/mnt/d/Documents/NGS/out/Manila_Bay/transients
# directory of kraken database
kraken_db=/mnt/d/Documents/NGS/kraken_db

# activate ngs conda environment (python version 3.8.13)
# all packages (except multiqc) were installed here: fastqc, bbtools, fastq_screen
conda activate ngs

# 1. check sequence quality of raw reads using fastQC
# create directory for fastQC results of raw reads
mkdir $out/fastqc/raw
# check quality of sequences using fastQC 
ls $in/*fastq.gz | parallel -j 4 "fastqc -o $out/fastqc/raw {}"

# 2. summarize fastQC results with multiQC
# change conda environment to use multiQC which runs with python 3.7
conda activate py3.7
# change directory to location of fastqc files
cd /mnt/d/Documents/NGS/out/Manila_Bay/fastqc
# run multiQC (https://multiqc.info/docs/#using-multiqc-reports)
# all 68 reads (34 samples and blanks)
multiqc $out/fastqc/raw/*_fastqc.zip -n /mnt/d/Documents/NGS/out/Manila_Bay/multiqc/fastqc_raw

# observations
# 1) per base sequence quality (sequence quality histograms)
# mostly R2 fail in mean quality score
# worse scores at the start of the read than at the end
# 2) per sequence quality scores
# all reads pass
# 3) per base sequence content
# biased sequence composition at the start of the reads
# balance starts at around base 15

# FASTQ generation already involves demultiplexing: https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html


# 3. prepare input files for merging / trimming
# change environment again for bbtools
conda activate ngs

# change directory to input (fastq) files
cd $in
# create list of Read 1 fastq
ls *_R1_* | sort > $transients/fastqR1
# create list of Read 2 fastq
ls *_R2_* | sort > $transients/fastqR2
# create tab-delimited file with full paths; this will be the input for bbmegere
paste $transients/fastqR1 $transients/fastqR2 > $transients/fastqIn

# Nextera adapters are not in the provided $bbmapResource/adapaters.fa 
# thus, we create our own reference with adapters and primers used for project
# check list of adapters written from nextera_adapters.R
cat $transients/addtl_nextera_adapters.fa
# contents of custom reference addtl_nextera_adapters.fa:
# 1) first round PCR primer overhangs 
# 2) locus-specific primers with overhangs
# 3) second round PCR primers (Nextera-style index primers where  i5 and i7 indicate the location of the barcode index sequences)
# i5 and i7 sequences added were from PGCs sequencing report
# P5-PCR index primer: 5’ AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC
# P7-PCR index primer: 5’ CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG

# 4. merge reads
# go to directory of input files
cd $in
# test parallel with input
parallel -j 4 --colsep '\t' echo "{2} {1}" :::: $transients/fastqIn
# merge reads 1 and 2
time parallel -j 1 --colsep '\t' --keep-order \
	bbmerge.sh "in1='{1}' in2='{2}'" \
	"out='$out/bbmerge/{=1 s/R1_001.fastq.gz/merged/;=}'" \
	"outu1='$out/bbmerge/{=1 s/R1_001.fastq.gz/R1_unmerged/;=}'" \
	"outu2='$out/bbmerge/{=2 s/R2_001.fastq.gz/R2_unmerged/;=}'" \
	"ordered=t outinsert='$out/bbmerge/bbmerge_insert_sizes'" \
	"outadapter='$out/bbmerge/bbmerge_consensus_adapter'" \
	"2>&1 | tee -a $out/bbmerge/log_bbmerge" :::: $transients/fastqIn

# create directory for fastQC results of merged/unmerged reads
mkdir $out/fastqc/merging
# run fastqc on merged sequences and unmerged sequences
ls $out/bbmerge/*_merged $out/bbmerge/*_unmerged | parallel -j 4 "fastqc -o $out/fastqc/merging {}"
# run multiqc 
conda activate py3.7
# for merged sequences:
multiqc $out/fastqc/merging/*_merged_fastqc.zip -n $out/multiqc/fastqc_merged
# for unmerged sequences:
multiqc $out/fastqc/merging/*_unmerged_fastqc.zip -n $out/multiqc/fastqc_unmerged

# 5. trimming of merged reads:
conda activate ngs
# go to directory of files to be trimmed
cd $out/bbmerge
# create list of merged reads to be trimmed
ls *_merged > $transients/fastqmerged
# run trimming of merged reads
time parallel -j 1 --keep-order \
	bbduk.sh "t=4 in='{1}'" \
	"out='$out/bbduk/{=1 s/merged/merged_trimmed/;=}'" \
	"minlength=80 ktrim=r k=23 mink=11 hdist=1 tbo tpe" \
	"ref='$transients/addtl_nextera_adapters.fa'" \
	"maxns=1 qtrim=r trimq=10" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_trim_fastqmerged" :::: $transients/fastqmerged
# notes:
# average no. of merged reads after trimming (maxn=1) is 25,621.47 (1459 - 94442 reads)
# difference in # of reads after trimming (between maxn=0 and maxn=1) ranges from 80 to 3045; average # of reads added is 727.7941 (0.95-12.95% increase, average 3.79%)
# hence, we make use of maxn=1

# 6. trimming of unmerged reads:
# create list of unmerged reads to be trimmed
sed 's/R1_001.fastq.gz/R1_unmerged/;s/R2_001.fastq.gz/R2_unmerged/' $transients/fastqIn > $transients/fastqunmerged
# run trimming of unmerged reads
time parallel -j 1 --colsep '\t' --keep-order \
	bbduk.sh "t=4 in1='{1}' in2='{2}'" \
	"out1='$out/bbduk/{=1 s/R1_unmerged/R1_unmerged_trimmed/;=}" \
	"out2='$out/bbduk/{=2 s/R2_unmerged/R2_unmerged_trimmed/;=}" \
	"minlength=80 ktrim=r k=23 mink=11 hdist=1 tbo tpe" \
	"ref='$transients/addtl_nextera_adapters.fa'" \
	"maxns=1 qtrim=r trimq=10" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_trim_fastqunmerged" :::: $transients/fastqunmerged

# create directory for fastQC results of merged/unmerged + trimmed reads
mkdir $out/fastqc/merging_trimming
# run fastqc on merged and unmerged reads that have been trimmed 
ls $out/bbduk/*_merged_trimmed $out/bbduk/*_unmerged_trimmed | parallel -j 4 "fastqc -o $out/fastqc/merging_trimming {}"
# run multiqc 
conda activate py3.7
# for merged, trimmed sequences:
multiqc $out/fastqc/merging_trimming/*_merged_trimmed_fastqc.zip -n $out/multiqc/fastqc_merged_trimmed
# for unmerged, trimmed sequences:
multiqc $out/fastqc/merging_trimming/*_unmerged_trimmed_fastqc.zip -n $out/multiqc/fastqc_unmerged_trimmed

# 7. trimming of raw reads (these will not undergo merging for comparison)
conda activate ngs
# go to directory of input files
cd $in
# list of reads to be trimmed is: $transients/fastqIn
time parallel -j 1 --colsep '\t' --keep-order \
	bbduk.sh "t=4 in1='{1}' in2='{2}'" \
	"out1='$out/bbduk/{=1 s/R1_001.fastq.gz/R1_trimmed/;=}" \
	"out2='$out/bbduk/{=2 s/R2_001.fastq.gz/R2_trimmed/;=}" \
	"minlength=80 ktrim=r k=23 mink=11 hdist=1 tbo tpe" \
	"ref='$transients/addtl_nextera_adapters.fa'" \
	"maxns=1 qtrim=r trimq=10" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_trim_fastqIn" :::: $transients/fastqIn
	
# create directory for fastQC results of trimmed reads
mkdir $out/fastqc/trimming
# run fastqc on trimmed reads
ls $out/bbduk/*[1-2]_trimmed | parallel -j 4 "fastqc -o $out/fastqc/trimming {}"
# run multiqc for trimmed sequences
conda activate py3.7
multiqc $out/fastqc/trimming/*_trimmed_fastqc.zip -n $out/multiqc/fastqc_trimmed

# 8. run fastq_screen to map reads against known genomes
# obtain reference genomes
fastq_screen --get_genomes --outdir $out/fastqscreen
# run fastq_screen for all trimmed sequences (merged and unmerged)
conda activate ngs
time fastq_screen --aligner bowtie2 --outdir "$out/fastqscreen" $(ls $out/bbduk/*_merged_trimmed $out/bbduk/*_unmerged_trimmed $out/bbduk/*[1-2]_trimmed) --conf /mnt/d/Documents/NGS/out/Manila_Bay/fastqscreen/FastQ_Screen_Genomes/fastq_screen.conf

# real    502m50.324s
# user    18m17.188s
# sys     78m39.906s

# multiqc of fastq_screen outputs
conda activate py3.7
# run multiqc for outputs of fastqscreen (merged, trimmed)
multiqc $out/fastqscreen/*_merged_trimmed_screen* -n $out/multiqc/fastqscreen_merged_trimmed
# run multiqc for outputs of fastqscreen (unmerged, trimmed)
multiqc $out/fastqscreen/*_unmerged_trimmed_screen* -n $out/multiqc/fastqscreen_unmerged_trimmed
# run multiqc for outputs of fastqscreen (trimmed only)
multiqc $out/fastqscreen/*[1-2]_trimmed_screen* -n $out/multiqc/fastqscreen_trimmed

# 9. filter reads
conda activate ngs
# go to directory of files to be filtered
cd $out/bbduk

# create list of merged + trimmed reads to filter
sed 's/$/_trimmed/' $transients/fastqmerged > $transients/fastqmerged_trimmed
# filter out phiX and sequencing artifacts from merged, trimmed reads
parallel -j 1 --keep-order \
	bbduk.sh "t=4 in='{1}' k=31" \
	"out='$out/bbduk/{=1 s/_merged_trimmed/_merged_trimmed_filtered/;=}'" \
	"ref='$bbmapResource/phix174_ill.ref.fa.gz,$bbmapResource/sequencing_artifacts.fa.gz'" \
	"stats='$out/bbduk/stats_filter_fastqmerged_trimmed' statscolumns=5" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_filter_fastqmerged_trimmed" :::: $transients/fastqmerged_trimmed
# for merged: only 1 result, 1B-16S-Metazoa_S23_L001_merged_trimmed_filtered

# create list of unmerged + trimmed reads to filter
sed 's/R1_unmerged/R1_unmerged_trimmed/;s/R2_unmerged/R2_unmerged_trimmed/' $transients/fastqunmerged > $transients/fastqunmerged_trimmed
# filter out phiX and sequencing artifacts from unmerged, trimmed reads
parallel -j 1 --colsep '\t' --keep-order \
	bbduk.sh "t=4 in1='{1}' in2='{2}' k=31" \
	"out1='$out/bbduk/{=1 s/R1_unmerged_trimmed/R1_unmerged_trimmed_filtered/;=}" \
	"out2='$out/bbduk/{=2 s/R2_unmerged_trimmed/R2_unmerged_trimmed_filtered/;=}" \
	"ref='$bbmapResource/phix174_ill.ref.fa.gz,$bbmapResource/sequencing_artifacts.fa.gz'" \
	"stats='$out/bbduk/stats_filter_fastqunmerged_trimmed' statscolumns=5" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_filter_fastqunmerged_trimmed" :::: $transients/fastqunmerged_trimmed
# result for unmerged reads: 
# 1B-16S-Metazoa_S23_L001_R1_unmerged_trimmed_filtered
# 1B-16S-Metazoa_S23_L001_R2_unmerged_trimmed_filtered

# create list of trimmed reads (did not undergo merging) to filter
sed 's/R1_001.fastq.gz/R1_trimmed/;s/R2_001.fastq.gz/R2_trimmed/' $transients/fastqIn > $transients/fastqtrimmed
# filter out phiX and sequencing artifacts from trimmed reads (no merging)
parallel -j 1 --colsep '\t' --keep-order \
	bbduk.sh "t=4 in1='{1}' in2='{2}' k=31" \
	"out1='$out/bbduk/{=1 s/R1_trimmed/R1_trimmed_filtered/;=}" \
	"out2='$out/bbduk/{=2 s/R2_trimmed/R2_trimmed_filtered/;=}" \
	"ref='$bbmapResource/phix174_ill.ref.fa.gz,$bbmapResource/sequencing_artifacts.fa.gz'" \
	"stats='$out/bbduk/stats_filter_fastqtrimmed' statscolumns=5" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_filter_fastqtrimmed" :::: $transients/fastqtrimmed
# filtered:
# 1B-16S-Metazoa_S23_L001_R1_trimmed_filtered  
# 1B-16S-Metazoa_S23_L001_R2_trimmed_filtered





# 10. coi kraken classification
# download taxonomy into $kraken_db directory (do only once)
kraken2-build --download-taxonomy --db $kraken_db

# add to library fasta files to be used
kraken2-build --add-to-library /mnt/d/Documents/NGS/entrez/kraken_fastas/worms_mares_coi.fasta --db $kraken_db --threads 4
# add note about contents of library so taxonomy directory can be re-used another time
echo 'Library contains /mnt/d/Documents/NGS/entrez/kraken_fastas/worms_mares_coi.fasta' > $kraken_db/readme

# build the kraken database
time kraken2-build --build --db $kraken_db --threads 4

# go to directory of input files for kraken
cd $out/bbduk

# create list of fastq for kraken (COI) classification (merged)
ls *UM*_merged_trimmed* > $transients/kraken_coi_merged
# kraken COI classification for merged (and trimmed) sequences
parallel -j 1 --keep-order kraken2 "--db $kraken_db --threads 4" \
	"--use-names --output $out/kraken/{=1 s/$/.kraken/;=} {1}" \
	"2>&1 | tee -a $out/kraken/log_kraken_coi_merged" :::: $transients/kraken_coi_merged
# create report
parallel -j 1 --keep-order kraken2 "--db $kraken_db --threads 4 --output -" \
	"--report $out/kraken/{=1 s/$/.krakenReport/;=} {1}" \
	"2>&1 | tee -a $out/kraken/log_krakenReport_coi_merged" :::: $transients/kraken_coi_merged
# use-names does not work for kraken reports

# create list of fastq for kraken (COI) classification (unmerged)	
# per read
cd $out/bbduk
ls *UM*R1_unmerged_trimmed | sort > $transients/kraken_coi_unmerged_R1
ls *UM*R2_unmerged_trimmed | sort > $transients/kraken_coi_unmerged_R2
# create tab-delimited file 
paste $transients/kraken_coi_unmerged_R1 $transients/kraken_coi_unmerged_R2 > $transients/kraken_coi_unmerged
# kraken COI classification for unmerged but paired sequences (trimmed)
parallel -j 1 --colsep '\t' --keep-order \
	kraken2 "--db $kraken_db --threads 4 --paired" \
	"--use-names --output $out/kraken/{=1 s/_R1_unmerged_trimmed/_unmerged_trimmed.kraken/;=} {1} {2}" \
	"2>&1 | tee -a $out/kraken/log_kraken_coi_unmerged" :::: $transients/kraken_coi_unmerged
# create report
parallel -j 1 --colsep '\t' --keep-order \
	kraken2 "--db $kraken_db --threads 4 --paired --output -" \
	"--report $out/kraken/{=1 s/_R1_unmerged_trimmed/_unmerged_trimmed.krakenReport/;=} {1} {2}" \
	"2>&1 | tee -a $out/kraken/log_krakenReport_coi_unmerged" :::: $transients/kraken_coi_unmerged

# create list of fastq for kraken (COI) classification (trimmed only)	
ls *UM*R1_trimmed | sort > $transients/kraken_coi_trimmed_R1
ls *UM*R2_trimmed | sort > $transients/kraken_coi_trimmed_R2
# create tab-delimited file 
paste $transients/kraken_coi_trimmed_R1 $transients/kraken_coi_trimmed_R2 > $transients/kraken_coi_trimmed
# kraken COI classification for trimmed reads (did not undergo merging)
parallel -j 1 --colsep '\t' --keep-order \
	kraken2 "--db $kraken_db --threads 4 --paired" \
	"--use-names --output $out/kraken/{=1 s/_R1_trimmed/_trimmed.kraken/;=} {1} {2}" \
	"2>&1 | tee -a $out/kraken/log_kraken_coi_trimmed" :::: $transients/kraken_coi_trimmed
# create report
parallel -j 1 --colsep '\t' --keep-order \
	kraken2 "--db $kraken_db --threads 4 --paired --output -" \
	"--report $out/kraken/{=1 s/_R1_trimmed/_trimmed.krakenReport/;=} {1} {2}" \
	"2>&1 | tee -a $out/kraken/log_krakenReport_coi_trimmed" :::: $transients/kraken_coi_trimmed

# 11. 16S kraken classification
# reuse taxonomy in $kraken_db: delete all files and directories except taxonomy
# within taxonomy directory, remove prelim_map

# add to library fasta files to be used
kraken2-build --add-to-library /mnt/d/Documents/NGS/entrez/kraken_fastas/worms_mares_16S_deduplicated.fasta --db $kraken_db --threads 4
# add note about contents of library so taxonomy directory can be re-used another time
echo 'Library contains /mnt/d/Documents/NGS/entrez/kraken_fastas/worms_mares_16S_deduplicated.fasta' > $kraken_db/readme

# build the kraken database
time kraken2-build --build --db $kraken_db --threads 4

# go to directory of input files for kraken
cd $out/bbduk

# create list of fastq for kraken (16S) classification (merged)
ls *16S*_merged_trimmed *OPH*_merged_trimmed > $transients/kraken_16S_merged
# kraken 16S classification for merged (and trimmed) sequences
parallel -j 1 --keep-order kraken2 "--db $kraken_db --threads 4" \
	"--use-names --output $out/kraken/{=1 s/$/.kraken/;=} {1}" \
	"2>&1 | tee -a $out/kraken/log_kraken_16S_merged" :::: $transients/kraken_16S_merged
# create report
parallel -j 1 --keep-order kraken2 "--db $kraken_db --threads 4 --output -" \
	"--report $out/kraken/{=1 s/$/.krakenReport/;=} {1}" \
	"2>&1 | tee -a $out/kraken/log_krakenReport_16S_merged" :::: $transients/kraken_16S_merged

# create list of fastq for kraken (16S) classification (unmerged)	
ls *16S*R1_unmerged_trimmed *OPH*R1_unmerged_trimmed | sort > $transients/kraken_16S_unmerged_R1
ls *16S*R2_unmerged_trimmed *OPH*R2_unmerged_trimmed | sort > $transients/kraken_16S_unmerged_R2
# create tab-delimited file 
paste $transients/kraken_16S_unmerged_R1 $transients/kraken_16S_unmerged_R2 > $transients/kraken_16S_unmerged
# kraken 16S classification for unmerged but paired sequences (trimmed)
parallel -j 1 --colsep '\t' --keep-order \
	kraken2 "--db $kraken_db --threads 4 --paired" \
	"--use-names --output $out/kraken/{=1 s/_R1_unmerged_trimmed/_unmerged_trimmed.kraken/;=} {1} {2}" \
	"2>&1 | tee -a $out/kraken/log_kraken_16S_unmerged" :::: $transients/kraken_16S_unmerged
# create report
parallel -j 1 --colsep '\t' --keep-order \
	kraken2 "--db $kraken_db --threads 4 --paired --output -" \
	"--report $out/kraken/{=1 s/_R1_unmerged_trimmed/_unmerged_trimmed.krakenReport/;=} {1} {2}" \
	"2>&1 | tee -a $out/kraken/log_krakenReport_16S_unmerged" :::: $transients/kraken_16S_unmerged

# create list of fastq for kraken (16S) classification (trimmed only)	
ls *16S*R1_trimmed *OPH*R1_trimmed | sort > $transients/kraken_16S_trimmed_R1
ls *16S*R2_trimmed *OPH*R2_trimmed | sort > $transients/kraken_16S_trimmed_R2
# create tab-delimited file 
paste $transients/kraken_16S_trimmed_R1 $transients/kraken_16S_trimmed_R2 > $transients/kraken_16S_trimmed
# kraken 16S classification for trimmed reads (did not undergo merging)
parallel -j 1 --colsep '\t' --keep-order \
	kraken2 "--db $kraken_db --threads 4 --paired" \
	"--use-names --output $out/kraken/{=1 s/_R1_trimmed/_trimmed.kraken/;=} {1} {2}" \
	"2>&1 | tee -a $out/kraken/log_kraken_16S_trimmed" :::: $transients/kraken_16S_trimmed
# create report
parallel -j 1 --colsep '\t' --keep-order \
	kraken2 "--db $kraken_db --threads 4 --paired --output -" \
	"--report $out/kraken/{=1 s/_R1_trimmed/_trimmed.krakenReport/;=} {1} {2}" \
	"2>&1 | tee -a $out/kraken/log_krakenReport_16S_trimmed" :::: $transients/kraken_16S_trimmed

# examine output using Pavian in R
pavian::runApp(port=5000)





# references
# https://uqbioinfo.github.io/pdf/2018-08-08-James.pdf 
# http://barcwiki.wi.mit.edu/wiki/SOPs/qc_shortReads
# https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/data-preprocessing/
# https://www.protocols.io/view/illumina-fastq-filtering-dm6gp9d5vzpn/v1?step=1

# counts per stage (summarized in read_counts.xlsx)
cd /mnt/d/Documents/NGS/input/fastqFiles/Manila_Bay
# get sample names
wc -l *R1*fastq.gz | awk '{print $2}' | cut -d '_' -f 1,2

# number of lines R1 raw fastq.gz (same as # of lines in R2)
# number of reads = number of lines / 4
# for-loop and zcat needed since files are compressed
for file in $(ls *R1*fastq.gz); do zcat $file | wc -l; done

# number of lines after merging
cd /mnt/d/Documents/NGS/out/Manila_Bay/bbmerge
wc -l *_merged | awk '{print $1}' 
# number of lines unmerged (# reads in R1 = # reads in R2)
wc -l *R1_unmerged | awk '{print $1}' 
# note: merged + unmerged = raw

# after merging + trimming
cd /mnt/d/Documents/NGS/out/Manila_Bay/bbduk 
# merged_trimmed
wc -l *_merged_trimmed | awk '{print $1}' 
# unmerged_trimmed
wc -l *R1_unmerged_trimmed | awk '{print $1}' 

# after merging + trimming + filtering
# merged_trimmed_filtered
wc -l *_merged_trimmed_filtered 
# unmerged_trimmed_filtered
wc -l *_unmerged_trimmed_filtered 

# after merging + trimming + classification
# number of lines here are equivalent to the number of reads
cd /mnt/d/Documents/NGS/out/Manila_Bay/kraken 
# merged_trimmed_classified
for file in $(ls *_merged_trimmed.kraken); do cut -f 1 $file | grep 'C' | wc -l; done
# unmerged_trimmed_classified
for file in $(ls *_unmerged_trimmed.kraken); do cut -f 1 $file | grep 'C' | wc -l; done

# counts processing without merging
cd /mnt/d/Documents/NGS/out/Manila_Bay/bbduk 
# trimmed
wc -l *_R1_trimmed | awk '{print $1}' 
# trimmed_filtered
wc -l *[1-2]_trimmed_filtered 

# after trimming + classification
# number of lines here are equivalent to the number of reads
cd /mnt/d/Documents/NGS/out/Manila_Bay/kraken 
# trimmed_classified
for file in $(ls *L001_trimmed.kraken); do cut -f 1 $file | grep 'C' | wc -l; done