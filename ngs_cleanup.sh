# set directory to where fastq files
cd /mnt/d/Documents/NGS/input/fastqFiles/Manila_Bay

# create directory for outputs
mkdir /mnt/d/Documents/NGS/out/Manila_Bay
# create directory for fastqc outputs
mkdir /mnt/d/Documents/NGS/out/Manila_Bay/fastqc
# create directory for multiqc outputs
mkdir /mnt/d/Documents/NGS/out/Manila_Bay/multiqc

# activate ngs conda environment (python version 3.8.13)
conda activate ngs

# check quality of sequences using fastQC 
ls *fastq.gz | parallel -j 4 "fastqc -o /mnt/d/Documents/NGS/out/Manila_Bay/fastqc {}"

# change conda environment to use multiQC
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
# 

# FASTQ generation already involves demultiplexing: https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html


