# set directory to location of fastq files
cd /mnt/d/Documents/NGS/input/fastqFiles/Manila_Bay

# create directory for outputs
mkdir /mnt/d/Documents/NGS/out/Manila_Bay
# create directory for fastqc outputs
mkdir /mnt/d/Documents/NGS/out/Manila_Bay/fastqc
# create directory for multiqc outputs
mkdir /mnt/d/Documents/NGS/out/Manila_Bay/multiqc
# create directory of transient files
mkdir /mnt/d/Documents/NGS/out/Manila_Bay/transients
# create directory for bbduk outputs
mkdir /mnt/d/Documents/NGS/out/Manila_Bay/bbduk
# create directory for fastq_screen outputs
mkdir /mnt/d/Documents/NGS/out/Manila_Bay/fastqscreen
# create directory for kraken outputs
mkdir /mnt/d/Documents/NGS/out/Manila_Bay/kraken

# activate ngs conda environment (python version 3.8.13)
conda activate ngs

# check quality of sequences using fastQC 
ls *fastq.gz | parallel -j 4 "fastqc -o /mnt/d/Documents/NGS/out/Manila_Bay/fastqc {}"

# change conda environment to use multiQC which runs with python 3.7
conda activate py3.7

# change directory to location of fastqc files
cd /mnt/d/Documents/NGS/out/Manila_Bay/fastqc
# run multiQC (https://multiqc.info/docs/#using-multiqc-reports)
# 1) all 68 reads (34 samples)
multiqc ./*_fastqc.zip -n /mnt/d/Documents/NGS/out/Manila_Bay/multiqc/fastqc_all
# 2) samples, uniminibarcode primer (36 reads or 18 samples)
multiqc ./[0-9]*UM*_fastqc.zip -n /mnt/d/Documents/NGS/out/Manila_Bay/multiqc/fastqc_samples_UM
# 3) samples, 16S primer (18 reads or 9 samples)
multiqc ./[0-9]*16S*_fastqc.zip -n /mnt/d/Documents/NGS/out/Manila_Bay/multiqc/fastqc_samples_16S
# 4) samples, OPH primer (6 reads or 3 samples
multiqc ./[0-9]*OPH*_fastqc.zip -n /mnt/d/Documents/NGS/out/Manila_Bay/multiqc/fastqc_samples_OPH
# 5) blanks + extraction negative (8 reads or 4 samples)
multiqc ./b*_fastqc.zip ./EN*_fastqc.zip -n /mnt/d/Documents/NGS/out/Manila_Bay/multiqc/fastqc_blanks
# 6) all samples (sans blanks and negative)
multiqc ./[0-9]*_fastqc.zip -n /mnt/d/Documents/NGS/out/Manila_Bay/multiqc/fastqc_samples_all

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





# change environment again for bbduk
conda activate ngs

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

# change directory to input files
cd $in
# create list of Read 1 fastq
ls *_R1_* | sort > $transients/fastqR1
# create list of Read 2 fastq
ls *_R2_* | sort > $transients/fastqR2
# create tab-delimited file with full paths
paste $transients/fastqR1 $transients/fastqR2 > $transients/fastqIn

# Nextera adapters are not in adapaters.fa reference so we append adapters and primers used for project
# check list of adapters written from nextera_adapters.R
cat $transients/addtl_nextera_adapters.fa
# append to adapters.fa to create custom adapters list
cat $bbmapResource/adapters.fa $transients/addtl_nextera_adapters.fa > $bbmapResource/custom_adapters.fa
# edit: will use addtl_nextera_adapters.fa later when trimming reads

# merge reads
# go to directory of input files
cd $in
# test parallel with input
parallel -j 4 --colsep '\t' echo "{2} {1}" :::: $transients/fastqIn
# merge reads 1 and 2
time parallel -j 1 --colsep '\t' --keep-order \
	bbmerge.sh "in1='{1}' in2='{2}'" \
	"out='$out/bbmerge/{=1 s/R1_001.fastq.gz/merged/;=}'" \
	"outu1='$out/bbmerge/{=1 s/R1_001.fastq.gz/unmerged_R1/;=}'" \
	"outu2='$out/bbmerge/{=2 s/R2_001.fastq.gz/unmerged_R2/;=}'" \
	"ordered=t outinsert='$out/bbmerge/bbmerge_insert_sizes'" \
	"outadapter='$out/bbmerge/bbmerge_consensus_adapter'" \
	"2>&1 | tee -a $out/bbmerge/log_bbmerge" :::: $transients/fastqIn

# run fastqc on merged sequences and unmerged sequences
ls $out/bbmerge/*_merged $out/bbmerge/*_unmerged_* | parallel -j 4 "fastqc -o $out/fastqc {}"
# run multiqc 
conda activate py3.7
# for merged sequences:
multiqc $out/fastqc/*merged_fastqc.zip -n $out/multiqc/fastqc_all_merged
# for unmerged sequences:
multiqc $out/fastqc/*_unmerged_*fastqc.zip -n $out/multiqc/fastqc_all_unmerged

# trim reads
conda activate ngs
# go to directory of files to be trimmed
cd $out/bbmerge

# trimming of merged reads:
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

# trimming of unmerged reads:
# create list of unmerged reads to be trimmed
sed 's/R1_001.fastq.gz/unmerged_R1/;s/R2_001.fastq.gz/unmerged_R2/' $transients/fastqIn > $transients/fastqunmerged
# run trimming of unmerged reads
time parallel -j 1 --colsep '\t' --keep-order \
	bbduk.sh "t=4 in1='{1}' in2='{2}'" \
	"out1='$out/bbduk/{=1 s/unmerged_R1/unmerged_trimmed_R1/;=}" \
	"out2='$out/bbduk/{=2 s/unmerged_R2/unmerged_trimmed_R2/;=}" \
	"minlength=80 ktrim=r k=23 mink=11 hdist=1 tbo tpe" \
	"ref='$transients/addtl_nextera_adapters.fa'" \
	"maxns=1 qtrim=r trimq=10" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_trim_fastqunmerged" :::: $transients/fastqunmerged

# run fastqc on merged, trimmed sequences and unmerged, trimmed sequences
ls $out/bbduk/*_merged_trimmed $out/bbduk/*_unmerged_trimmed_R* | parallel -j 4 "fastqc -o $out/fastqc {}"
# run multiqc 
conda activate py3.7
# for merged, trimmed sequences:
multiqc $out/fastqc/*merged_trimmed_fastqc.zip -n $out/multiqc/fastqc_all_merged_trimmed
# for unmerged, trimmed sequences:
multiqc $out/fastqc/*_unmerged_trimmed_R*fastqc.zip -n $out/multiqc/fastqc_all_unmerged_trimmed

# fastq_screen to map reads against known genomes
# obtain reference genomes
fastq_screen --get_genomes --outdir $out/fastqscreen
# run fastq_screen for all trimmed sequences (merged and unmerged)
conda activate ngs
time fastq_screen --aligner bowtie2 --outdir "$out/fastqscreen" $(ls $out/bbduk/*merged_trimmed $out/bbduk/*_unmerged_trimmed*) --conf /mnt/d/Documents/NGS/out/Manila_Bay/fastqscreen/FastQ_Screen_Genomes/fastq_screen.conf
# real    315m57.701s
# user    11m28.859s
# sys     52m15.016s

# multiqc of fastq_screen outputs
conda activate py3.7
# run multiqc for outputs of fastqscreen (merged, trimmed)
multiqc $out/fastqscreen/*merged_trimmed_screen* -n $out/multiqc/fastqscreen_merged_trimmed
# run multiqc for outputs of fastqscreen (unmerged, trimmed)
multiqc $out/fastqscreen/*_unmerged_trimmed_R*_screen* -n $out/multiqc/fastqscreen_unmerged_trimmed

# filter reads
conda activate ngs
# create list of merged trimmed reads to filter
sed 's/$/_trimmed/' $transients/fastqmerged > $transients/fastqmerged_trimmed
# go to directory of files to be filtered
cd $out/bbduk
# filter out phiX and sequencing artifacts from merged, trimmed reads
parallel -j 1 --keep-order \
	bbduk.sh "t=4 in='{1}' k=31" \
	"out='$out/bbduk/{=1 s/merged_trimmed/merged_trimmed_filtered/;=}'" \
	"ref='$bbmapResource/phix174_ill.ref.fa.gz,$bbmapResource/sequencing_artifacts.fa.gz'" \
	"stats='$out/bbduk/stats_filter_fastqmerged_trimmed' statscolumns=5" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_filter_fastqmerged_trimmed" :::: $transients/fastqmerged_trimmed
# initial run using unmerged reads: only 1B-16S-Metazoa_S23_L001_bbdukOut2_R1  1B-16S-Metazoa_S23_L001_bbdukOut2_R2 with filtered when paired reads were input
# for merged: only 1 result, 1B-16S-Metazoa_S23_L001_merged_trimmed_filtered

# create list of unmerged trimmed reads to filter
sed 's/unmerged_R1/unmerged_trimmed_R1/;s/unmerged_R2/unmerged_trimmed_R2/' $transients/fastqunmerged > $transients/fastqunmerged_trimmed
# filter out phiX and sequencing artifacts from unmerged, trimmed reads
parallel -j 1 --colsep '\t' --keep-order \
	bbduk.sh "t=4 in1='{1}' in2='{2}' k=31" \
	"out1='$out/bbduk/{=1 s/unmerged_trimmed_R1/unmerged_trimmed_filtered_R1/;=}" \
	"out2='$out/bbduk/{=2 s/unmerged_trimmed_R2/unmerged_trimmed_filtered_R2/;=}" \
	"ref='$bbmapResource/phix174_ill.ref.fa.gz,$bbmapResource/sequencing_artifacts.fa.gz'" \
	"stats='$out/bbduk/stats_filter_fastqunmerged_trimmed' statscolumns=5" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_filter_fastqunmerged_trimmed" :::: $transients/fastqunmerged_trimmed
# result for unmerged reads: 1B-16S-Metazoa_S23_L001_unmerged_trimmed_filtered_R1  1B-16S-Metazoa_S23_L001_unmerged_trimmed_filtered_R2







# kraken classification
# download taxonomy into $kraken_db directory
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
# kraken COI classification for merged sequences
parallel -j 1 --keep-order kraken2 "--db $kraken_db --threads 4" \
	"--use-names --output $out/kraken/{=1 s/$/.kraken/;=} {1}" \
	"2>&1 | tee -a $out/kraken/log_kraken_coi_merged" :::: $transients/kraken_coi_merged
# create report
parallel -j 1 --keep-order kraken2 "--db $kraken_db --threads 4 --output -" \
	"--report $out/kraken/{=1 s/$/.krakenReport/;=} {1}" \
	"2>&1 | tee -a $out/kraken/log_krakenReport_coi_merged" :::: $transients/kraken_coi_merged
# use-names does not work for kraken reports

# create list of fastq for kraken (COI) classification (unmerged but paired)	
# per read
cd $out/bbduk
ls *UM*unmerged_trimmed_R1 | sort > $transients/kraken_coi_unmerged_R1
ls *UM*unmerged_trimmed_R2 | sort > $transients/kraken_coi_unmerged_R2
# create tab-delimited file 
paste $transients/kraken_coi_unmerged_R1 $transients/kraken_coi_unmerged_R2 > $transients/kraken_coi_unmerged
# kraken COI classification for unmerged but paired sequences
parallel -j 1 --colsep '\t' --keep-order \
	kraken2 "--db $kraken_db --threads 4 --paired" \
	"--use-names --output $out/kraken/{=1 s/_R1/.kraken/;=} {1} {2}" \
	"2>&1 | tee -a $out/kraken/log_kraken_coi_unmerged" :::: $transients/kraken_coi_unmerged
# create report
parallel -j 1 --colsep '\t' --keep-order \
	kraken2 "--db $kraken_db --threads 4 --paired --output -" \
	"--report $out/kraken/{=1 s/_R1/.krakenReport/;=} {1} {2}" \
	"2>&1 | tee -a $out/kraken/log_krakenReport_coi_unmerged" :::: $transients/kraken_coi_unmerged




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
wc -l *-UM*R1* | awk '{print $2}' | cut -d '_' -f 1,2
wc -l *-16S*R1* | awk '{print $2}' | cut -d '_' -f 1,2
wc -l *-OPH*R1* | awk '{print $2}' | cut -d '_' -f 1,2

# R1 raw files 
for file in $(ls *-UM*R1*); do zcat $file | wc -l; done
for file in $(ls *-16S*R1*); do zcat $file | wc -l; done
for file in $(ls *-OPH*R1*); do zcat $file | wc -l; done
# R2 raw files (just the same as R1)
for file in $(ls *-UM*R2*); do zcat $file | wc -l; done
for file in $(ls *-16S*R2*); do zcat $file | wc -l; done
for file in $(ls *-OPH*R2*); do zcat $file | wc -l; done


# merged
cd /mnt/d/Documents/NGS/out/Manila_Bay/bbmerge
wc -l *-UM*merged | awk '{print $1}' 
wc -l *-16S*merged | awk '{print $1}' 
wc -l *-OPH*merged | awk '{print $1}' 

# merged_trimmed
cd /mnt/d/Documents/NGS/out/Manila_Bay/bbduk 
wc -l *-UM*merged_trimmed | awk '{print $1}' 
wc -l *-16S*merged_trimmed | awk '{print $1}' 
wc -l *-OPH*merged_trimmed | awk '{print $1}' 

# merged_trimmed_filtered
wc -l *-UM*merged_trimmed_filtered | awk '{print $1}' 
wc -l *-16S*merged_trimmed_filtered | awk '{print $1}' 
wc -l *-OPH*merged_trimmed_filtered | awk '{print $1}' 

# unmerged
cd /mnt/d/Documents/NGS/out/Manila_Bay/bbmerge
wc -l *-UM*unmerged_R1 | awk '{print $1}' 
wc -l *-16S*unmerged_R1 | awk '{print $1}' 
wc -l *-OPH*unmerged_R1 | awk '{print $1}' 
wc -l *-UM*unmerged_R2 | awk '{print $1}' 
wc -l *-16S*unmerged_R2 | awk '{print $1}' 
wc -l *-OPH*unmerged_R2 | awk '{print $1}' 

# unmerged_trimmed
cd /mnt/d/Documents/NGS/out/Manila_Bay/bbduk 
wc -l *-UM*unmerged_trimmed_R1| awk '{print $1}' 
wc -l *-16S*unmerged_trimmed_R1 | awk '{print $1}' 
wc -l *-OPH*unmerged_trimmed_R1 | awk '{print $1}' 
wc -l *-UM*unmerged_trimmed_R2| awk '{print $1}' 
wc -l *-16S*unmerged_trimmed_R2 | awk '{print $1}' 
wc -l *-OPH*unmerged_trimmed_R2 | awk '{print $1}' 

# unmerged_trimmed_filtered
wc -l *-UM*unmerged_trimmed_filtered_R1 | awk '{print $1}' 
wc -l *-16S*unmerged_trimmed_filtered_R1 | awk '{print $1}' 
wc -l *-OPH*unmerged_trimmed_filtered_R1 | awk '{print $1}' 
wc -l *-UM*unmerged_trimmed_filtered_R2 | awk '{print $1}' 
wc -l *-16S*unmerged_trimmed_filtered_R2 | awk '{print $1}' 
wc -l *-OPH*unmerged_trimmed_filtered_R2 | awk '{print $1}' 

# merged_trimmed_classified
cd /mnt/d/Documents/NGS/out/Manila_Bay/kraken 
for file in $(ls *-UM*_merged_trimmed.kraken); do cut -f 1 $file | grep 'C' | wc -l; done
# unmerged_trimmed_classified
for file in $(ls *-UM*_unmerged_trimmed.kraken); do cut -f 1 $file | grep 'C' | wc -l; done

