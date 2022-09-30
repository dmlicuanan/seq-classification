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
export out=/mnt/d/Documents/NGS/out/Manila_Bay
# in = directory of raw reads
export in=/mnt/d/Documents/NGS/input/fastqFiles/Manila_Bay
# directory of bbmapresources
export bbmapResource=/home/dmlicuanan/miniconda3/envs/ngs/opt/bbmap-38.18/resources
# directory of transient files
export transients=/mnt/d/Documents/NGS/out/Manila_Bay/transients

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

# merge reads
# go to directory of input files
cd $in
# test parallel with input
parallel -j 4 --colsep '\t' echo "{2} {1}" :::: $transients/fastqIn
# merge reads 1 and 2
parallel -j 4 --colsep '\t' --keep-order \
	bbmerge.sh "t=4 in1='{1}' in2='{2}'" \
	"out='$out/bbmerge/{=1 s/R1_001.fastq.gz/merged/;=}'" \
	"outu1='$out/bbmerge/{=1 s/R1_001.fastq.gz/unmerged_R1/;=}'" \
	"outu2='$out/bbmerge/{=2 s/R2_001.fastq.gz/unmerged_R2/;=}'" \
	"outinsert='$out/bbmerge/bbmerge_insert_sizes'" \
	"outadapter='$out/bbmerge/bbmerge_consensus_adapter'" \
	"2>&1 | tee -a $out/bbmerge/log_bbmerge" :::: $transients/fastqIn

# run fastqc on merged sequences
ls $out/bbmerge/*_merged | parallel -j 4 "fastqc -o $out/fastqc {}"
# run multiqc for merged sequences
conda activate py3.7
multiqc $out/fastqc/*merged_fastqc.zip -n $out/multiqc/fastqc_all_merged

# trim reads
conda activate ngs
# go to directory of files to be trimmed
cd $out/bbmerge
# create list of merged reads to be trimmed
ls *_merged > $transients/fastqmerged
# run trimming 
time parallel -j 4 --keep-order \
	bbduk.sh "t=4 in='{1}'" \
	"out='$out/bbduk/{=1 s/merged/merged_trimmed/;=}'" \
	"minlength=80 ktrim=r k=23 mink=11 hdist=1 tbo tpe" \
	"ref='$transients/addtl_nextera_adapters.fa'" \
	"maxns=1 qtrim=r trimq=10" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_trim" :::: $transients/fastqmerged

# fastq_screen to map reads against known genomes
# obtain reference genomes
fastq_screen --get_genomes --outdir $out/fastqscreen
# run fastq_screen
fastq_screen --aligner bowtie2 --outdir "$out/fastqscreen" $(ls $out/bbduk/*merged_trimmed) --conf /mnt/d/Documents/NGS/out/Manila_Bay/fastqscreen/FastQ_Screen_Genomes/fastq_screen.conf

# multiqc of fastq_screen outputs
cd $out/fastqscreen
# activate py3.7
conda activate py3.7
# run multiqc for outputs of fastqscreen
multiqc ./*merged_trimmed_screen* -n $out/multiqc/fastqscreen_merged_trimmed

# filter reads
conda activate ngs
# create list of merged trimmed reads to filter
sed 's/$/_trimmed/' $transients/fastqmerged > $transients/fastqmerged_trimmed
# go to directory of merged trimmed reads to be filtered
cd $out/bbduk
# filter out phiX and sequencing artifacts
parallel -j 4 --keep-order \
	bbduk.sh "t=4 in='{1}' k=31" \
	"out='$out/bbduk/{=1 s/merged_trimmed/merged_trimmed_filtered/;=}'" \
	"ref='$bbmapResource/phix174_ill.ref.fa.gz,$bbmapResource/sequencing_artifacts.fa.gz'" \
	"stats='$out/bbduk/stats_filter' statscolumns=5" \
	"2>&1 | tee -a $out/bbduk/log_bbduk_filter" :::: $transients/fastqmerged_trimmed

# only 1B-16S-Metazoa_S23_L001_bbdukOut2_R1  1B-16S-Metazoa_S23_L001_bbdukOut2_R2 with filtered


# references
# https://uqbioinfo.github.io/pdf/2018-08-08-James.pdf 
# http://barcwiki.wi.mit.edu/wiki/SOPs/qc_shortReads
# https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/data-preprocessing/
# https://www.protocols.io/view/illumina-fastq-filtering-dm6gp9d5vzpn/v1?step=1
