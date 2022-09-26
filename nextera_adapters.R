# script for generating fasta containing adapters
# sequences will be used for trimming reads

# read in library validation and indexing report table from sequencing report
lib <- read.table("D:/Documents/Repositories/seq_classification/data/library_indexing.txt", col.names = c("sample_no", "sample_name", "conc", "lib_size", "N7XX", "i7", "S5XX", "i5"))

# sequences will be based on https://dnatech.genomecenter.ucdavis.edu/faqs/how-to-prepare-samples-for-multiplexed-amplicon-sequencing-on-illumina-systems/
# 1) first round PCR primer overhangs 
# Forward overhang P5-tag: 5’ TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG-[locus-specific sequence]
# Reverse overhang P7-tag: 5’ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-[locus-specific sequence]
overhangs <- data.frame(name = c("Nextera_forward_overhang", "Nextera_reverse_overhang"), primer = c("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"))

# 2) second round PCR primers (Nextera-style index primers where  i5 and i7 indicate the location of the barcode index sequences)
# P5-PCR index primer: 5’ AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC
# P7-PCR index primer: 5’ CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG
# according to: https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2019/03/illumina-adapter-sequences-2019-1000000002694-10.pdf:
# N7XX = index 1 (i7) adapters
# S5XX = index 2 (i5) adapters

# unique values of N7XX
n7 <- unique(lib[,c("N7XX", "i7")])
# paste other parts of primer 
n7$primer <- paste0("CAAGCAGAAGACGGCATACGAGAT", n7$i7, "GTCTCGTGGGCTCGG")
# add name
n7$name <- paste0("Nextera_P7_PCR_index_primer_", n7$N7XX)

# unique values of S5XX
s5 <- unique(lib[,c("S5XX", "i5")])
# paste other parts of primer 
s5$primer <- paste0("AATGATACGGCGACCACCGAGATCTACAC", s5$i5, "TCGTCGGCAGCGTC")
# add name
s5$name <- paste0("Nextera_P5_PCR_index_primer_", s5$S5XX)

# combine tables
all <- rbind(overhangs, n7[,c(4,3)], s5[,c(4,3)])

# load writeFasta() function
# write .fas file from a vector of IDs (accession numbers info etc.) and sequences
writeFasta <- function(id, seq, file){
  # make a dummy output
  writeTemp <- 1:(2*length(id))
  # replace odd numbers with length
  writeTemp[seq(1, length(writeTemp), 2)] <- paste(">", id, sep = "")
  # replace even numbers with sequences
  writeTemp[seq(2, length(writeTemp), 2)] <- seq
  # write the dummy as fasta
  write(writeTemp, file = file)
  return(NULL)
}

# write fasta that will be appended to adapters.fa in bbmap resources directory
writeFasta(id = all$name, seq = all$primer, 
           file = "D:/Documents/NGS/out/Manila_Bay/transients/addtl_nextera_adapters.fa")