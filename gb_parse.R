# script for parsing GenBank flatfiles
# contains functions for creating fasta files from GenBank flatfiles
# functions should work for files that contain multiple sequences
# functions be able to extract the COI portion if mitochondrial genomes are provided

# GB flatfile information: 
# https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#Top

# working directory: path of .gb files
setwd("D:/Documents/NGS/entrez/wormstaxlist_gb/compiled_gb")
# check number of sequences expected
# grep "^LOCUS" *_ncbi.gb -w | wc -l = 470806 > 470915 
# grep "^DEFINITION" *_ncbi.gb -w | wc -l = 470806 > 470915
# grep "^ORIGIN" *_ncbi.gb -w | wc -l = 470794 > 470904
# grep "^ORIGIN\|^WGS \|"^CONTIG" *_ncbi.gb -w | wc -l = 470915
# grep "^//" *_ncbi.gb -w | wc -l = 470803 > 470915

# view possible formats
readLines("Aaptos_aaptos_ncbi.gb") # 1 COI sequence
readLines("Aaptos_papillata_ncbi.gb") # 2 COI sequences
readLines("Anemonia_viridis_ncbi.gb", 705) # complete genome

# list of .gbs to process = 31708
files <- list.files(pattern = "_ncbi.gb")

# because of inconsistencies in number of sequences expected, need to check which files have patterns outside expectations
check <- function(gbfile) {
  # read file
  gb <- readLines(gbfile)
  
  # record indices of markers which should occur once in each record
  m1 <- grep("^LOCUS", gb)
  m2 <- grep("^DEFINITION", gb)
  m3 <- grep("^ORIGIN", gb)
  m4 <- grep("^//", gb)
  # compile in list
  occ <- list(m1, m2, m3, m4)
  lens <- sapply(occ, length)
  
  # return file name if length of occurrences are not equal
  if (min(lens) == max(lens)) { NULL } 
  # else { return(gbfile) }
  else { return(lens) } # for bad files
}
# apply function over entire list of files (453.63 seconds)
badfiles <- unlist(lapply(files, check))
# apply function over bad files to check occurrences per marker
# edit else part of check function
badmat <- matrix(unlist(lapply(badfiles, check)), byrow = TRUE, ncol = 4)
colSums(badmat)
# consistent with grep counts in Ubuntu:
# ORIGIN lacking by 12 from LOCUS AND DEFINITION
# // lacking by 3

# create list of files to re-download from NCBI
bf <- gsub("_ncbi.gb", " [ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])", badfiles)
bf <- gsub("_", " ", bf)
# write file readable in Unix
# f <- file("D:/Documents/NGS/entrez/wormstaxlist_badfiles", open="wb")
# cat(bf, file = f, sep = "\n")
# close(f)

ptm <- proc.time()
badfiles <- unlist(lapply(files, check))
proc.time() - ptm
beep(2)


temp <- vector()
for (i in 1:length(files)) {
  gb <- readLines(files[i])
  temp[i] <- length(gb)
  gbtofas(gb)
}
  
for (i in 1:length(files)) {
  gb <- readLines(files[i])
  temp[i] <- length(gb)
  gbtofas(gb)
}
  




gb <- readLines("Aaptos_aaptos_ncbi.gb")
genbankFile <- gb
gettaxaInfo(gb)
getseqInfo(gb)


getseqInfo <- function(genbankFile){
  temp <- grep(x = genbankFile, pattern = "ORIGIN|//")
  seqInfo <- genbankFile[(temp[1] + 1):(temp[2] - 1)]
  rm(temp)
  seqInfo <- paste(gsub(x = seqInfo, pattern = "[[:digit:]]| ",
                        replacement = ""),
                   collapse = "")
  return(seqInfo)
}

substr()