# script for parsing GenBank flatfiles
# contains functions for creating fasta files from GenBank flatfiles
# functions should work for files that contain multiple sequences
# functions be able to extract the COI portion if mitochondrial genomes are provided

# GB flatfile information: 
# https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#Top

# working directory: path of .gb files
setwd("D:/Documents/NGS/entrez/wormstaxlist_gb/compiled_gb")
# check number of sequences expected
# count after > is number after re-downloading bad files
# grep "^LOCUS" *_ncbi.gb -w | wc -l = 470806 > 470915 
# grep "^DEFINITION" *_ncbi.gb -w | wc -l = 470806 > 470915
# grep "^ORIGIN" *_ncbi.gb -w | wc -l = 470794 > 470904
# grep "^ORIGIN\|^WGS \|"^CONTIG" *_ncbi.gb -w | wc -l = 470915
# grep "^//" *_ncbi.gb -w | wc -l = 470803 > 470915
# grep "^     source" *_ncbi.gb -w | wc -l = 470915
# grep "db_xref=\"taxon:" -n *_ncbi.gb | wc -l = 470915

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
  
  # record indices of markers
  m1 <- grep("^LOCUS", gb)
  m2 <- grep("^DEFINITION", gb)
  m3 <- grep("^ORIGIN", gb)
  m4 <- grep("^//", gb)
  m5 <- grep("^WGS |^CONTIG", gb) # replace ORIGIN in records with no sequences
  # compile in list
  occ <- list(m1, m2, m3, m4)
  lens <- sapply(occ, length)
  
  # return file name and counts per marker if length of occurrences are not equal
  if (min(lens) == max(lens)) { NULL } 
  else { return(c(gbfile, lens, length(m5))) } 
}

# apply function over entire list of files (453.63 seconds)
badfiles <- unlist(lapply(files, check))

# split vector
badmat <- as.data.frame(matrix(badfiles, byrow = TRUE, ncol = 6))
# format dataframe
colnames(badmat) <- c("file", "locus", "def", "ori", "end", "wgs")
badmat[,-1] <- lapply(badmat[,-1], as.numeric)
# check counts
all(badmat$locus == badmat$def)
all(badmat$locus == badmat$end)
all(badmat$locus != badmat$ori)
all(badmat$locus == (badmat$ori + badmat$wgs))
colSums(badmat[,-1])
# consistent with grep counts in Ubuntu:
# ORIGIN is not a consistent marker in each record
# WGS or CONTIG replace ORIGIN, but there are no sequences  

# create list of files to re-download from NCBI
# bf <- gsub("_ncbi.gb", " [ORGN] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])", badmat$file)
# bf <- gsub("_", " ", bf)
# write file readable in Unix
# f <- file("D:/Documents/NGS/entrez/wormstaxlist_badfiles", open="wb")
# cat(bf, file = f, sep = "\n")
# close(f)

# 1) function for converting .gb (gbfile) to .fasta (fastafile)
gbtofasta <- function(gbfile, fastafile) {
  # read file
  gb <- readLines(gbfile)
  
  # record start and end markers per record
  start <- grep("^LOCUS", gb)
  end <- grep("^//", gb)
  origin <- grep("^ORIGIN", gb)
  # empty vector
  f <- vector()
  
  for (i in 1:length(start)) {
    # is there an ORIGIN marker in between LOCUS and //?
    if (any(origin %in% start[i]:end[i])) { 
      # restrict actions to record
      record <- gb[start[i]:end[i]]
      
      # get accession version
      version <- record[grep("^VERSION", record)]
      version <- trimws(gsub("^VERSION", "", version))
      
      # get definition
      def <- record[grep("^DEFINITION", record):(grep("^ACCESSION", record)-1)]
      def <- paste(trimws(gsub("^DEFINITION", "", def)), collapse = " ")
      # remove last period
      def <- sub("[.]$", "", def)
      
      # get sequence data
      seq <- record[(grep("^ORIGIN", record)+1):(grep("^//", record)-1)]
      seq <- toupper(paste(gsub("[[:digit:]]| ", "", seq), collapse = ""))
      
      # assemble fasta
      id <- paste0(">", version, " ", def)
      f <- c(f, id, seq)
    }
    else { NULL }
  }
  # write fasta if vector not empty
  if (length(f) > 0) { write(f, fastafile) }
  else { NULL }
}

# file for testing gbtofasta function 
# gbfile <- "Aspergillus_fumigatus_ncbi.gb"

# convert all .gb files to .fasta (69.16483 elapsed mins)
for (i in 1:length(files)) {
  # vector of fasta file names
  f <- sub("_ncbi.gb", "_ncbi.fasta", files)
  # paste path with file name
  f <- paste0("D:/Documents/NGS/entrez/wormstaxlist_gb/wormtaxlist_fasta", "/", f)
  # apply gbtofasta function
  gbtofasta(files[i], f[i])
}







# test file 
gbfile <- "Aspergillus_fumigatus_ncbi.gb"
gbfile <- 
i <- 1 # no ORIGIN
i <- 2
i <- 3
i <- 4 # COI

# 2) function for converting .gb (gbfile) to .fasta (fastafile)
# this time selecting sequence that corresponds to COI
gbtofasta <- function(gbfile, fastafile) {
  # read file
  gb <- readLines(gbfile)
  
  # record start and end markers per record
  start <- grep("^LOCUS", gb)
  end <- grep("^//", gb)
  origin <- grep("^ORIGIN", gb)
  # empty vector
  f <- vector()
  
  for (i in 1:length(start)) {
    # is there an ORIGIN marker in between LOCUS and //?
    if (any(origin %in% start[i]:end[i])) { 
      # restrict actions to record
      record <- gb[start[i]:end[i]]
      
      # get accession version
      version <- record[grep("^VERSION", record)]
      version <- trimws(gsub("^VERSION", "", version))
      
      # get definition
      def <- record[grep("^DEFINITION", record):(grep("^ACCESSION", record)-1)]
      def <- paste(trimws(gsub("^DEFINITION", "", def)), collapse = " ")
      # remove last period
      def <- sub("[.]$", "", def)
      
      # restrict to info under FEATURES (include line of ORIGIN)
      ind <- grep("^FEATURES|^ORIGIN", record)
      feats <- record[(ind[1]+1):(ind[2])]
      
      # get info under source
      # indices of sections; feats as reference 
      indSec <- grep("^     [[:alpha:]]|^ORIGIN", feats)
      # index of section with "source"; section indices as reference
      indSource <- grep("^     source", feats[indSec])
      # contents of source section
      sub <- feats[indSec[indSource]:(indSec[indSource+1]-1)]
      # get taxid
      taxid <- gsub("\\\"|\\s|/", "", sub[grep("db_xref=\"taxon:", sub)])
      taxid <- sub("db_xref=taxon:", "", taxid)
      # get location of whole sequence
      sourceLoc <- 
      
      # get sections with gene; section indices as reference
      indGene <- grep("^     gene", feats[indSec])
      
      # examines lines per gene feature
      for (j in 1:length(indGene)) {
        sub <- feats[indSec[indGene][j]:(indSec[indGene+1]-1)[j]]
        pattern <- "/gene.*CO1|/gene.*COI|/gene.*COX1|/gene.*COXI|/gene.*CO-I|/gene.*CO-1"

        if (any(grepl(pattern, sub, ignore.case = TRUE))) { 
          loc <- gsub("gene|\\s", "", sub[grep("^     gene", sub)])
          gene <- gsub("/gene=|\\s|\"", "", sub[grep("/gene=", sub)])}
        else { NULL }
        
      }
      
      
      indSec[indGene]
      indSec[indGene+1]-1
      feats[indSec[indGene]]
        
      feats[grep("^                     /", feats)]
      
      # get sequence data
      # seq <- record[(grep("^ORIGIN", record)+1):(grep("^//", record)-1)]
      # seq <- toupper(paste(gsub("[[:digit:]]| ", "", seq), collapse = ""))
      
      # assemble fasta
      # id <- paste0(">", version, " ", def)
      # f <- c(f, id, seq)
    }
    else { NULL }
  }
  # write fasta if vector not empty
  if (length(f) > 0) { write(f, fastafile) }
  else { NULL }
}




ptm <- proc.time()
# put function here:

proc.time() - ptm
beep(5)



