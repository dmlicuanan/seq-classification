# script for parsing GenBank flatfiles
# contains functions for creating fasta files from GenBank flatfiles
# also contains functions that parse files with multiple .gb records, extract sequence information, and extract the COI portion of sequences if mitogenomes are provided

# references:
# https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#Top
# https://www.insdc.org/submitting-standards/feature-table/#1

# CONVERSION OF COI GENBANK FILES TO FASTAS
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

# function for converting .gb (gbfile) to .fasta (fastafile)
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


# some downloaded sequences do not contain COI exclusively
# need to cut out COI sequences only
# function descriptions:
# 1) gbParse: takes .gb file and separates records
# 2) recInfo: takes record, returns list of definition, version, taxid, source (range of sequence)
# 3) coiInfo: takes record, returns gene name and location

# 1) gbParse
# input .gb file which contains many sequences
# output list of records
gbParse <- function(gbfile) {
  # read .gb file
  gb <- readLines(gbfile)
  # mark start and end markers per record
  start <- grep("^LOCUS", gb)
  end <- grep("^//", gb)
  # marks the presence of a sequence
  origin <- grep("^ORIGIN", gb)
  
  # empty list
  recs <- list()
  
  # only proceed if no. of start markers and no. of end markers is equal
  if (length(start) == length(end)) {
    for (i in 1:length(start)) {
      # is there an ORIGIN marker in between LOCUS and //?
      # if yes, proceed, since the record holds a sequence
      if (any(origin %in% start[i]:end[i])) {
        # hold in ith element of list
        recs[[i]] <- gb[start[i]:end[i]]
      } else { print(paste("gbParse: there is no sequence data in record #", i, "of", gbfile)) }
    }
  } else { print(paste("gbParse: number of start and end markers not equal in", gbfile)) }
  # return list sans NULL elements
  return(recs[lengths(recs) != 0])
}

# 2) recInfo
# function to extract the following from each record:
# definition, version, taxid, source (range of sequence)
recInfo <- function(gbrecord) {
  # create empty list where to store record information
  list <- vector("list", length = 4)
  names(list) <- c("def", "version", "taxid", "source")
  
  # get accession version
  list$version <- gbrecord[grep("^VERSION", gbrecord)]
  list$version <- trimws(gsub("^VERSION", "", list$version))
  
  # get definition
  # line index containing start of definition
  ind1 <- grep("^DEFINITION", gbrecord)
  # indices of headers
  headers <- grep("^[[:alpha:]]", gbrecord)
  # line index of end of definition
  ind2 <- headers[which(headers == ind1) + 1] - 1
  list$def <- gbrecord[ind1:ind2]
  list$def <- paste(trimws(gsub("^DEFINITION", "", list$def)), collapse = " ")
  # remove last period
  list$def <- sub("[.]$", "", list$def)
  
  # line index containing "source"
  ind1 <- grep("^     source", gbrecord)
  # indices of feature starts
  headers <- grep("^     [[:alpha:]]", gbrecord)
  # line index of end of source key
  ind2 <- headers[which(headers %in% ind1) + 1] - 1
  
  # is there only one source key in the record?
  if (length(ind1) == 1) {
    # restrict actions to section under source key
    sub <- gbrecord[ind1:ind2]
    # get taxid
    list$taxid <- gsub("\\\"|\\s|/", "", sub[grep("db_xref=\"taxon:", sub)])
    list$taxid <- sub("db_xref=taxon:", "", list$taxid)
    # replace taxid with NA, and print message
    if (length(list$taxid) == 0) { 
      list$taxid <- NA
      print(paste("recInfo: sequence tagged", list$version, "has no taxid"))}
    
    # get location range of sequence
    list$source <- gsub("source|\\s", "", gbrecord[ind1])
    
  } else { print(paste("recInfo: there is more than one source key in sequence tagged", list$version))}
  
  # output list of record information 
  return(list)
}

# 3) coiInfo
# function to extract gene and location information of COI
coiInfo <- function(gbrecord) {
  # create empty vectors where locations and COI gene names will be compiled
  loc <- vector()
  gene <- vector()
  
  # get features with "gene"
  ind1 <- grep("^     gene", gbrecord)
  # line indices where all features start + index of ORIGIN 
  headers <- grep("^     [[:alpha:]]|^ORIGIN", gbrecord)
  # indices where each "gene" feature ends
  ind2 <- headers[which(headers %in% ind1) + 1] - 1
  
  for (j in 1:length(ind1)) {
    # examines lines per gene feature
    sub <- gbrecord[ind1[j]:ind2[j]]
    # COI pattern to search for 
    pattern <- "/gene.*CO1\"|/gene.*COI\"|/gene.*COX1\"|/gene.*COXI\"|/gene.*CO-I\"|/gene.*CO-1\""
    
    # if any COI pattern is in the gene feature, extract location and gene name
    if (any(grepl(pattern, sub, ignore.case = TRUE))) { 
      # which lines contain the start of location
      stLine <- grep("^     gene", sub)
      # which lines contain the end of location
      endLine <- grep("^                     /", sub)[1] - 1
      # location line
      l <- sub[stLine:endLine]
      # remove large spaces and collapse as one line
      l <- gsub("^                     ", "", l)
      l <- paste(l, collapse = "")
      
      # get location
      loc <- c(loc, gsub("gene|\\s", "", l))
      # get gene
      gene <- c(gene, gsub("/gene=|\\s|\"", "", sub[grep("/gene=", sub)]))
      }
  }
  # if no COI pattern found in sections, print message in console
  if (length(gene) == 0) { 
    # get accession version
    ver <- gbrecord[grep("^VERSION", gbrecord)]
    ver <- trimws(gsub("^VERSION", "", ver))
    # print message
    print(paste("geneInfo: no COI gene information in sequence tagged", ver)) }
  
  # put gene and loc in dataframe to check duplicates
  d <- unique(data.frame(loc, gene))
  # re-assign back to vectors
  loc <- d$loc
  gene <- d$gene
  # return named list
  return(list(gene = gene, loc = loc))
}

# location data is not simply encoded, so need to check possible characters that arise in location info 
# all .gb files
files <- list.files(pattern = "_ncbi.gb")
# run gbParse, recInfo, coiInfo (32.44633 elapsed mins; 38.01567 elapsed mins)
for (i in 1:length(files)) {
  # parse .gb files
  gblist <- gbParse(files[i])
  # nested loop
  for (j in 1:length(gblist)) {
    # perform actions per record
    record <- gblist[[j]]
    # get record and gene info
    l1 <- recInfo(record)
    l2 <- coiInfo(record)
    
    # base nrows on number of genes listed as COI
    df <- data.frame(gene = l2$gene, loc = l2$loc)
    # add unique data per record
    df$def <- l1$def
    df$version <- l1$version
    df$taxid <- l1$taxid
    df$source <- l1$source
    df$file <- files[i]
    
    # append df to file
    # 1st run -- when coiInfo only took one line of loc info
    # write.table(df, "D:/Documents/NGS/entrez/wormstaxlist_gb/locInfo.txt", append = TRUE, col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(df, "D:/Documents/NGS/entrez/wormstaxlist_gb/locInfo_2.txt", append = TRUE, col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  }
}
# notes -- consistent with grep counts!
# 11 .gb records with no sequence
# [1] "gbParse: there is no sequence data in record # 1 of Aspergillus_fumigatus_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 35 of Clunio_marinus_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 1 of Dirofilaria_immitis_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 1 of Hortaea_werneckii_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 28 of Perca_fluviatilis_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 29 of Perca_fluviatilis_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 1 of Pichia_spartinae_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 2 of Pichia_spartinae_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 12 of Saccharomyces_cerevisiae_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 125 of Saccharomyces_cerevisiae_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 8 of Scatophagus_argus_ncbi.gb"

# check location information 
# locRaw <- read.table("D:/Documents/NGS/entrez/wormstaxlist_gb/locInfo.txt", sep = "\t", quote = "", comment.char = "") # 471011 rows
locRaw <- read.table("D:/Documents/NGS/entrez/wormstaxlist_gb/locInfo_2.txt", sep = "\t", quote = "", comment.char = "") # 471009 rows
# add names to columns
names(locRaw) <- c("gene", "loc", "def", "version", "taxid", "source", "file")

# check the format of source
# https://www.insdc.org/submitting-standards/feature-table/#1
# confirm if it always starts with 1..(length of sequence)
temp <- locRaw$source
# do all begin with 1..
all(grepl("^1\\.\\.", temp))
# do all follow the format: 1..digit
all(grepl("^1\\.\\.[[:digit:]]*$", temp))
unique(sub("^1\\.\\.[[:digit:]]*$", "", temp))
# sequence lengths
len <- as.integer(sub("^1\\.\\.", "", temp))
# range of sequence lengths: 20 to 678653
range(len)
# distribution of lengths
hist(len, breaks = 100)
hist(sort(len)[1:450000], breaks = 100)
# 94.56% of lengths are less than 1000 bp
sum(len < 1000) / length(len)

# which accession #s are repeated?
df <- data.frame(table(locRaw$version))
# vector of accession numbers that appear more than once
dups <- df[df$Freq > 1,]$Var1
# print duplicates based on accession #
temp <- locRaw[locRaw$version %in% dups,]
temp <- temp[order(temp$version),]
head(temp[, c("version", "taxid", "source", "file")])
# note: some duplicates are from species synonyms
# e.g.: Cynarina lacrymalis is synonymous with Sclerophyllia margariticola
# both are linked to the same .gb record AB117246.1

# number of unique accession #s in df showing duplicates
length(unique(temp$version)) # 3754
# number of unique entries based on record contents alone
nrow(unique(temp[,c('def', 'version', 'taxid', 'source')]))
# note: number of version is indicative of the #s of unique records

# remove duplicates due to species synonyms to see other duplicates 
# duplication will then be due to multiple COI locations per record
# remove file column
loc <- locRaw[, -7]
# remove duplicates based on all remaining variables
loc <- loc[!duplicated(loc), ]
# check again which accession numbers are repeated
df <- data.frame(table(loc$version))
# vector of accession numbers that appear more than once
dups <- as.character(df[df$Freq > 1,]$Var1)
# print duplicates based on accession #
temp <- loc[loc$version %in% dups,]
temp <- temp[order(temp$version),]
temp[, c("version", "source", "loc", "gene")]
# there are 97 records associated with more than one COI gene location 
length(unique(temp$version))

# coiChecker
# function for printing all coi features in .gb record if there is more than one
coiChecker <- function(gbrecord) {
  # create empty vectors where locations and COI gene names will be compiled
  loc <- vector()
  gene <- vector()
  # vector where feature section will be stored
  feature <- vector()
  
  # get features with "gene" or "CDS"
  ind1 <- grep("^     gene", gbrecord)
  # line indices where all features start + index of ORIGIN 
  headers <- grep("^     [[:alpha:]]|^ORIGIN", gbrecord)
  # indices where each "gene" feature ends
  ind2 <- headers[which(headers %in% ind1) + 1] - 1
  
  for (j in 1:length(ind1)) {
    # examines lines per gene feature
    sub <- gbrecord[ind1[j]:ind2[j]]
    # COI pattern to search for 
    pattern <- "/gene.*CO1\"|/gene.*COI\"|/gene.*COX1\"|/gene.*COXI\"|/gene.*CO-I\"|/gene.*CO-1\""
    
    # if any COI pattern is in the gene feature, extract location and gene name
    if (any(grepl(pattern, sub, ignore.case = TRUE))) { 
      loc <- c(loc, gsub("gene|\\s", "", sub[grep("^     gene", sub)]))
      gene <- c(gene, gsub("/gene=|\\s|\"", "", sub[grep("/gene=", sub)]))
      # save feature in vector
      feature <- c(feature, sub)
    }
  }
  
  # get accession version
  ver <- gbrecord[grep("^VERSION", gbrecord)]
  ver <- trimws(gsub("^VERSION", "", ver))
  
  # if no COI pattern found in sections, print message in console
  if (length(gene) == 0) { 
    # print message
    print(paste("geneInfo: no COI gene information in sequence tagged", ver)) }
  
  # if there is more than one COI feature, print accession # and feature
  if (length(gene) > 1) {
    print(ver)
    print(feature)
  }
}

# vector of accessions associated with more than one COI location
vers <- unique(temp$version)
# find file names linked to example accession #s
f <- unique(locRaw[locRaw$version %in% vers, ]$file) 
# mapping of accession version with file name
unique(locRaw[locRaw$version %in% vers, c("version", "file")])

# print sections when COI locations > 1
for (i in 1:length(f)) {
  gblist <- gbParse(f[i])
  
  for (j in 1:length(gblist)) {
    # perform actions per record
    record <- gblist[[j]]
    coiChecker(record)
    }
}
# [1] "MK562756.1"
# [1] "     gene            91095..92591"  
# [2] "                     /gene=\"cox1\""
# [3] "     gene            286304..287800"
# [4] "                     /gene=\"cox1\""
# [1] "FJ429092.1"
# [1] "     gene            2095..3627"   
# [2] "                     /gene=\"CO1\""
# [3] "     gene            8548..10080"  
# [4] "                     /gene=\"CO1\""
# [1] "KC701764.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8549..10081"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701763.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701762.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701761.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701760.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701759.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701757.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701755.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701754.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701753.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701752.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701751.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701750.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701749.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701748.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701747.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701746.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701745.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701744.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8549..10081"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701743.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701742.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701741.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701740.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701738.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701737.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701736.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701735.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701734.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701733.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701732.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701731.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701730.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701729.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701728.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701727.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701725.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8548..10080"   
# [4] "                     /gene=\"COX1\""
# [1] "KC701724.1"
# [1] "     gene            2095..3627"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8557..10089"   
# [4] "                     /gene=\"COX1\""
# [1] "NC_011581.1"
# [1] "     gene            2095..3627"                 
# [2] "                     /gene=\"COX1\""             
# [3] "                     /db_xref=\"GeneID:7042918\""
# [4] "     gene            8548..10080"                
# [5] "                     /gene=\"COX1\""             
# [6] "                     /db_xref=\"GeneID:7042913\""
# [1] "NC_016423.1"
# [1] "     gene            1416..2948"                  
# [2] "                     /gene=\"COX1\""              
# [3] "                     /db_xref=\"GeneID:11452066\""
# [4] "     gene            15055..16587"                
# [5] "                     /gene=\"COX1\""              
# [6] "                     /db_xref=\"GeneID:11452077\""
# [1] "AP012225.1"
# [1] "     gene            1416..2948"   
# [2] "                     /gene=\"CO1\""
# [3] "     gene            15055..16587" 
# [4] "                     /gene=\"CO1\""
# [1] "JQ062883.1"
# [1] "     gene            12447..19832"  
# [2] "                     /gene=\"cox1\""
# [3] "     gene            13440..15491"  
# [4] "                     /gene=\"cox1\""
# [1] "JQ062882.1"
# [1] "     gene            12445..19830"  
# [2] "                     /gene=\"cox1\""
# [3] "     gene            13438..15489"  
# [4] "                     /gene=\"cox1\""
# [1] "JQ062881.1"
# [1] "     gene            11838..19223"  
# [2] "                     /gene=\"cox1\""
# [3] "     gene            12081..14882"  
# [4] "                     /gene=\"cox1\""
# [5] "     gene            12831..14882"  
# [6] "                     /gene=\"cox1\""
# [1] "NC_040118.1"
# [1] "     gene            1433..2965"                  
# [2] "                     /gene=\"COX1\""              
# [3] "                     /locus_tag=\"EJ553_mgp17\""  
# [4] "                     /db_xref=\"GeneID:38574330\""
# [5] "     gene            9092..10624"                 
# [6] "                     /gene=\"COX1\""              
# [7] "                     /locus_tag=\"EJ553_mgp11\""  
# [8] "                     /db_xref=\"GeneID:38574321\""
# [1] "MG833837.1"
# [1] "     gene            1433..2965"    
# [2] "                     /gene=\"cox1\""
# [3] "     gene            9092..10624"   
# [4] "                     /gene=\"cox1\""
# [1] "KT809323.1"
# [1] "     gene            complement(295..373)"
# [2] "                     /gene=\"cox1\""      
# [3] "                     /pseudo"             
# [4] "     gene            13604..15169"        
# [5] "                     /gene=\"cox1\""      
# [1] "NC_016465.1"
# [1] "     gene            complement(80..1645)"        
# [2] "                     /gene=\"COX1\""              
# [3] "                     /db_xref=\"GeneID:11473170\""
# [4] "     gene            15119..16684"                
# [5] "                     /gene=\"COX1\""              
# [6] "                     /db_xref=\"GeneID:11473165\""
# [1] "JN700935.1"
# [1] "     gene            complement(80..1645)"
# [2] "                     /gene=\"cox1\""      
# [3] "     gene            15119..16684"        
# [4] "                     /gene=\"cox1\""      
# [1] "FO082257.1"
# [1] "     gene            24414..39193"                    
# [2] "                     /gene=\"cox1\""                  
# [3] "     gene            <24669..27116"                   
# [4] "                     /gene=\"cox1\""                  
# [5] "                     /locus_tag=\"cpurp_mito_Cox1-1\""
# [6] "     gene            <30561..31238"                   
# [7] "                     /gene=\"cox1\""                  
# [8] "                     /locus_tag=\"cpurp_mito_Cox1-7\""
# [9] "     gene            <32093..32893"                   
# [10] "                     /gene=\"cox1\""                  
# [11] "                     /locus_tag=\"cpurp_mito_Cox1-9\""
# [1] "NC_047242.1"
# [1] "     gene            complement(38111..53468)"    
# [2] "                     /gene=\"cox1\""              
# [3] "                     /locus_tag=\"HJF13_mgp060\"" 
# [4] "                     /db_xref=\"GeneID:54599960\""
# [5] "     gene            complement(53591..54520)"    
# [6] "                     /gene=\"cox1\""              
# [7] "                     /locus_tag=\"HJF13_mgp047\"" 
# [8] "                     /db_xref=\"GeneID:54599973\""
# [1] "MN978926.1"
# [1] "     gene            complement(38111..53468)"
# [2] "                     /gene=\"cox1\""          
# [3] "     gene            complement(53591..54520)"
# [4] "                     /gene=\"cox1\""          
# [1] "MW592987.1"
# [1] "     gene            8344..10408"   
# [2] "                     /gene=\"cox1\""
# [3] "     gene            8554..8821"    
# [4] "                     /gene=\"cox1\""
# [5] "     gene            9464..9703"    
# [6] "                     /gene=\"cox1\""
# [1] "MW592989.1"
# [1] "     gene            3335..5419"    
# [2] "                     /gene=\"cox1\""
# [3] "     gene            3542..3818"    
# [4] "                     /gene=\"cox1\""
# [5] "     gene            4461..4705"    
# [6] "                     /gene=\"cox1\""
# [1] "EU068697.1"
# [1] "     gene            1425..2957"   
# [2] "                     /gene=\"CO1\""
# [3] "     gene            8454..9986"   
# [4] "                     /gene=\"CO1\""
# [1] "NC_009734.1"
# [1] "     gene            1425..2957"                 
# [2] "                     /gene=\"COX1\""             
# [3] "                     /db_xref=\"GeneID:5469443\""
# [4] "     gene            8454..9986"                 
# [5] "                     /gene=\"COX1\""             
# [6] "                     /db_xref=\"GeneID:5469447\""
# [1] "MW450849.1"
# [1] "     gene            2057..3589"    
# [2] "                     /gene=\"cox1\""
# [3] "     gene            8454..9986"    
# [4] "                     /gene=\"cox1\""
# [1] "MW592986.1"
# [1] "     gene            1..3645"       
# [2] "                     /gene=\"cox1\""
# [3] "     gene            893..2974"     
# [4] "                     /gene=\"cox1\""
# [1] "LN901196.1"
# [1] "     gene            complement(9..87)"
# [2] "                     /gene=\"cox1\""   
# [3] "                     /pseudo"          
# [4] "     gene            13444..15009"     
# [5] "                     /gene=\"cox1\""   
# [1] "LN901197.1"
# [1] "     gene            complement(111..189)"
# [2] "                     /gene=\"cox1\""      
# [3] "                     /pseudo"             
# [4] "     gene            13552..15117"        
# [5] "                     /gene=\"cox1\""      
# [1] "NC_026908.1"
# [1] "     gene            1987..3519"                  
# [2] "                     /gene=\"COX1\""              
# [3] "                     /locus_tag=\"YC02_gp16\""    
# [4] "                     /db_xref=\"GeneID:24143643\""
# [5] "     gene            9016..10548"                 
# [6] "                     /gene=\"COX1\""              
# [7] "                     /locus_tag=\"YC02_gp10\""    
# [8] "                     /db_xref=\"GeneID:24143649\""
# [1] "KP336702.1"
# [1] "     gene            1987..3519"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            9016..10548"   
# [4] "                     /gene=\"COX1\""
# [1] "JN700945.1"
# [1] "     gene            complement(<1..798)"
# [2] "                     /gene=\"cox1\""     
# [3] "     gene            14213..15778"       
# [4] "                     /gene=\"cox1\""     
# [1] "NC_063457.1"
# [1] "     gene            1..1476"                     
# [2] "                     /gene=\"cox1\""              
# [3] "                     /locus_tag=\"NDC64_mgp01\""  
# [4] "                     /note=\"copy1\""             
# [5] "                     /db_xref=\"GeneID:72631909\""
# [6] "     gene            1495..3591"                  
# [7] "                     /gene=\"cox1\""              
# [8] "                     /locus_tag=\"NDC64_mgp02\""  
# [9] "                     /note=\"copy2\""             
# [10] "                     /db_xref=\"GeneID:72631842\""
# [1] "MW849268.1"
# [1] "     gene            1..1476"        
# [2] "                     /gene=\"cox1\"" 
# [3] "                     /note=\"copy1\""
# [4] "     gene            1495..3591"     
# [5] "                     /gene=\"cox1\"" 
# [6] "                     /note=\"copy2\""
# [1] "MW849269.1"
# [1] "     gene            1..2004"         
# [2] "                     /gene=\"cox1\""  
# [3] "                     /note=\"copy 1\""
# [4] "     gene            2816..4576"      
# [5] "                     /gene=\"cox1\""  
# [6] "                     /note=\"copy 2\""
# [7] "     gene            5145..9745"      
# [8] "                     /gene=\"cox1\""  
# [9] "                     /note=\"copy 3\""
# [1] "JN700948.1"
# [1] "     gene            complement(<1..948)"
# [2] "                     /gene=\"cox1\""     
# [3] "     gene            14348..>15194"      
# [4] "                     /gene=\"cox1\""     
# [1] "NC_020348.1"
# [1] "     gene            2041..3573"                  
# [2] "                     /gene=\"COX1\""              
# [3] "                     /gene_synonym=\"COI\""       
# [4] "                     /db_xref=\"GeneID:14658336\""
# [5] "     gene            15321..16853"                
# [6] "                     /gene=\"COX1\""              
# [7] "                     /gene_synonym=\"COI\""       
# [8] "                     /db_xref=\"GeneID:14658337\""
# [1] "AB715401.1"
# [1] "     gene            2041..3573"   
# [2] "                     /gene=\"COI\""
# [3] "     gene            15321..16853" 
# [4] "                     /gene=\"COI\""
# [1] "NC_031832.1"
# [1] "     gene            complement(16336..17946)"    
# [2] "                     /gene=\"cox1\""              
# [3] "                     /locus_tag=\"BOO99_gp058\""  
# [4] "                     /db_xref=\"GeneID:30214182\""
# [5] "     gene            58929..60539"                
# [6] "                     /gene=\"cox1\""              
# [7] "                     /locus_tag=\"BOO99_gp015\""  
# [8] "                     /db_xref=\"GeneID:30214195\""
# [1] "AP017433.1"
# [1] "     gene            complement(16336..17946)"
# [2] "                     /gene=\"cox1\""          
# [3] "     gene            58929..60539"            
# [4] "                     /gene=\"cox1\""          
# [1] "MW592988.1"
# [1] "     gene            10577..14861"  
# [2] "                     /gene=\"cox1\""
# [3] "     gene            10754..12963"  
# [4] "                     /gene=\"cox1\""
# [5] "     gene            14038..14564"  
# [6] "                     /gene=\"cox1\""
# [1] "KT809328.1"
# [1] "     gene            complement(9..85)"
# [2] "                     /gene=\"cox1\""   
# [3] "                     /pseudo"          
# [4] "     gene            <13268..14817"    
# [5] "                     /gene=\"cox1\""   
# [1] "LN901209.1"
# [1] "     gene            complement(53..129)"
# [2] "                     /gene=\"cox1\""     
# [3] "                     /pseudo"            
# [4] "     gene            13296..14861"       
# [5] "                     /gene=\"cox1\""     
# [1] "LN901210.1"
# [1] "     gene            complement(227..304)"
# [2] "                     /gene=\"cox1\""      
# [3] "                     /pseudo"             
# [4] "     gene            13638..15203"        
# [5] "                     /gene=\"cox1\""      
# [1] "gbParse: there is no sequence data in record # 12 of Saccharomyces_cerevisiae_ncbi.gb"
# [1] "gbParse: there is no sequence data in record # 125 of Saccharomyces_cerevisiae_ncbi.gb"
# [1] "HG994156.1"
# [1] "     gene            11349..21552"                                              
# [2] "                     /gene=\"COX1\""                                            
# [3] "                     /locus_tag=\"C2U11_6151\""                                 
# [4] "     gene            order(11349..11517,13966..14039,15555..16033,17044..17223,"
# [5] "                     18692..18897,20172..20195,21080..21552)"                   
# [6] "                     /gene=\"COX1\""                                            
# [7] "                     /locus_tag=\"C2U11_6147\""                                 
# [8] "     gene            order(11349..11517,13966..14039,15555..16020,16021..16033,"
# [9] "                     17044..17223,18692..18897,20172..20195,21080..21552)"      
# [10] "                     /gene=\"COX1\""                                            
# [11] "                     /locus_tag=\"C2U11_6143\""                                 
# [1] "LR999888.1"
# [1] "     gene            11349..21552"                                              
# [2] "                     /gene=\"COX1\""                                            
# [3] "                     /locus_tag=\"EO220_6159\""                                 
# [4] "     gene            order(11349..11517,13966..14039,15555..16033,17044..17223,"
# [5] "                     18692..18897,20172..20195,21080..21552)"                   
# [6] "                     /gene=\"COX1\""                                            
# [7] "                     /locus_tag=\"EO220_6151\""                                 
# [8] "     gene            order(11349..11517,13966..14039,15555..16033,17044..17223,"
# [9] "                     18692..18897,20172..20195,21080..21552)"                   
# [10] "                     /gene=\"COX1\""                                            
# [11] "                     /locus_tag=\"EO220_6155\""                                 
# [1] "MT471321.1"
# [1] "     gene            74680..76260"              
# [2] "                     /gene=\"cox1\""            
# [3] "     gene            complement(459681..461261)"
# [4] "                     /gene=\"cox1\""            
# [1] "MT661575.1"
# [1] "     gene            1421..2953"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            8445..9977"    
# [4] "                     /gene=\"COX1\""
# [1] "KT962062.1"
# [1] "     gene            1421..2953"           
# [2] "                     /gene=\"COX1\""       
# [3] "                     /gene_synonym=\"COI\""
# [4] "     gene            8445..9977"           
# [5] "                     /gene=\"COX1\""       
# [6] "                     /gene_synonym=\"COI\""
# [1] "EU658923.1"
# [1] "     gene            1423..2955"   
# [2] "                     /gene=\"CO1\""
# [3] "     gene            8445..9977"   
# [4] "                     /gene=\"CO1\""
# [1] "EU660577.1"
# [1] "     gene            1434..2966"   
# [2] "                     /gene=\"CO1\""
# [3] "     gene            8461..9993"   
# [4] "                     /gene=\"CO1\""
# [1] "EU660576.1"
# [1] "     gene            1421..2953"   
# [2] "                     /gene=\"CO1\""
# [3] "     gene            8443..9975"   
# [4] "                     /gene=\"CO1\""
# [1] "NC_010636.1"
# [1] "     gene            1421..2953"                 
# [2] "                     /gene=\"COX1\""             
# [3] "                     /db_xref=\"GeneID:6261943\""
# [4] "     gene            8443..9975"                 
# [5] "                     /gene=\"COX1\""             
# [6] "                     /db_xref=\"GeneID:6261929\""
# [1] "NC_058301.1"
# [1] "     gene            2034..3566"                  
# [2] "                     /gene=\"COX1\""              
# [3] "                     /locus_tag=\"LI390_mgp17\""  
# [4] "                     /db_xref=\"GeneID:68212244\""
# [5] "     gene            15404..16936"                
# [6] "                     /gene=\"COX1\""              
# [7] "                     /locus_tag=\"LI390_mgp05\""  
# [8] "                     /db_xref=\"GeneID:68212255\""
# [1] "MT733875.1"
# [1] "     gene            2034..3566"    
# [2] "                     /gene=\"cox1\""
# [3] "     gene            15404..16936"  
# [4] "                     /gene=\"cox1\""
# [1] "NC_006354.1"
# [1] "     gene            1400..2938"                 
# [2] "                     /gene=\"COX1\""             
# [3] "                     /db_xref=\"GeneID:3101905\""
# [4] "     gene            8423..9955"                 
# [5] "                     /gene=\"COX1\""             
# [6] "                     /db_xref=\"GeneID:3101913\""
# [1] "AB240153.1"
# [1] "     gene            2035..3567"   
# [2] "                     /gene=\"CO1\""
# [3] "     gene            15269..16801" 
# [4] "                     /gene=\"CO1\""
# [1] "AB158364.1"
# [1] "     gene            1400..2938"   
# [2] "                     /gene=\"COI\""
# [3] "     gene            8423..9955"   
# [4] "                     /gene=\"COI\""
# [1] "JX473257.1"
# [1] "     gene            <1..327"                                       
# [2] "                     /gene=\"coxI\""                                
# [3] "                     /note=\"cytochrome c oxidase subunit I; numt\""
# [4] "                     /pseudo"                                       
# [5] "     gene            572..>1093"                                    
# [6] "                     /gene=\"coxI\""                                
# [7] "                     /note=\"cytochrome c oxidase subunit I; numt\""
# [8] "                     /pseudo"                                       
# [1] "KT020766.1"
# [1] "     gene            complement(33..329)"
# [2] "                     /gene=\"cox1\""     
# [3] "                     /pseudo"            
# [4] "     gene            13636..15393"       
# [5] "                     /gene=\"COX1\""     
# [1] "NC_031213.1"
# [1] "     gene            complement(33..329)"         
# [2] "                     /gene=\"COX1\""              
# [3] "                     /locus_tag=\"BI101_gp01\""   
# [4] "                     /pseudo"                     
# [5] "                     /db_xref=\"GeneID:29061427\""
# [6] "     gene            13636..15393"                
# [7] "                     /gene=\"COX1\""              
# [8] "                     /locus_tag=\"BI101_gp02\""   
# [9] "                     /db_xref=\"GeneID:29061421\""
# [1] "KJ845633.1"
# [1] "     gene            2016..3548"    
# [2] "                     /gene=\"COX1\""
# [3] "     gene            15094..16626"  
# [4] "                     /gene=\"COX1\""
# [1] "NC_007893.1"
# [1] "     gene            2016..3548"                 
# [2] "                     /gene=\"COX1\""             
# [3] "                     /db_xref=\"GeneID:3974256\""
# [4] "     gene            15098..16630"               
# [5] "                     /gene=\"COX1\""             
# [6] "                     /db_xref=\"GeneID:3974257\""
# [1] "AB240152.1"
# [1] "     gene            2016..3548"   
# [2] "                     /gene=\"CO1\""
# [3] "     gene            15098..16630" 
# [4] "                     /gene=\"CO1\""
# [1] "AB086202.1"
# [1] "     gene            1384..2916"    
# [2] "                     /gene=\"cox1\""
# [3] "     gene            8395..9927"    
# [4] "                     /gene=\"cox1\""

# inspect sequence spans with location operators
d <- unique(locRaw$loc)
d[grep("[[:alpha:]]", d)] # 8616 / 471009 or 2%
# print sequence spans with  nested location operators
d[grep("\\(.*\\(", d)]
# [1] "complement(order(7564..8430,9532..10395))"                               
# [2] "complement(join(16980..17024,1..1548))"                                  
# [3] "join(3434..4153,5046..5056,complement(12235..12327),21866..21911,22081..22815)"
# [4] "complement(join(15808..15855,1..1542))"

# inspect sequence spans with no location operators
d <- unique(locRaw$loc)
d <- d[grep("[[:alpha:]]", d, invert = TRUE)]
# removing < and >, do all sequence spans follow the format of digit..digit
t <- gsub("<|>", "", d)
all(grepl("^[[:digit:]]*\\.\\.[[:digit:]]*$", t))
# range of sequence lengths
start <- as.integer(sub("\\.\\.[[:digit:]].*", "", t))
end <- as.integer(sub("[[:digit:]].*\\.\\.", "", t))
len <- end - start
hist(sort(len)[1:5700], breaks = 100)
range(len) 
# find sequences with length less than 100 bp
smallspans <- d[which(len < 100)]
locRaw[locRaw$loc %in% smallspans,-c(3,5)]
# some sequences are short since the record is an amalgamation of different regions.
# e.g. KJ666477.1 Goniastrea pectinata voucher MTQ G61908 cytochrome c oxidase subunit I (COI) gene, partial cds; COI-trnM intergenic spacer, complete sequence; and tRNA-Met (trnM) gene, partial sequence; mitochondrial. whose COI span is <1..46
# there are also sequences that are less than 10 bp:
locRaw[locRaw$loc %in% d[which(len < 10)],-c(3,5)]
# e.g.: LN824143.1 Anabathron hedleyi mitochondrial partial COI gene for cytochrome oxidase subunit 1, specimen voucher AMS, C.467225. with COI span <1..>3, but source span is 1..625, and no other genes are under features

# % of records whose COI location span is less than 100 = 0.305% = 1436/471009
(nrow(locRaw[locRaw$loc %in% smallspans,]) / nrow(locRaw))*100
# % of sequences with small spans whose location span == source span; 55/1436
t <- locRaw[locRaw$loc %in% smallspans,]
t$loc <- gsub("<|>", "", t$loc)
t[which(t$loc == t$source), -c(3,5)]
# % of sequences whose COI span == source span is 93.33%
t <- locRaw
t$loc <- gsub("<|>", "", t$loc)
(nrow(t[t$loc == t$source,]) / nrow(t))*100

# function for extracting sequence from GenBank record
# assumes ORIGIN marker exists in the record
seqData <- function(gbrecord) {
  seq <- gbrecord[(grep("^ORIGIN", gbrecord)+1):(grep("^//", gbrecord)-1)]
  seq <- toupper(paste(gsub("[[:digit:]]| ", "", seq), collapse = ""))
  return(seq)
}

# function for cutting a given sequence based on a location span
# seq = sequence to be cut
# start = vector of sequence starts
# end = vector of sequence stops
seqCut <- function(seq, start, end) {
  # vector to hold cut sequences
  vec <- vector()
  for (i in 1:length(start)) {
    vec[i] <- substr(seq, start = start[i], stop = end[i])
  }
  return(vec)
}

# files for conversion to COI fastas
files <- list.files(pattern = "_ncbi.gb")

# for loop for converting .gbs to .fastas
# load required package:
library(seqinr)
# sequences are restricted to COI locations
# 41.44767 elapsed minutes to run (NCBI format)
# 46.162 elapsed minutes to run (kraken format)
for (i in 1:length(files)) {
  gblist <- gbParse(files[i])
  
  # empty vector where fasta will be assembled per .gb
  f <- vector()
  # do actions per record in the .gb file
  for (j in 1:length(gblist)) {
    # take jth record
    record <- gblist[[j]]
    
    # get record metadata 
    data <- recInfo(record)
    
    # id info for sequence (unique for each record)
    # id <- paste0(">", data$version, " ", data$def)
    
    # id if NCBI format will be overriden for kraken format
    id <- paste0(">", "kraken:taxid|", data$taxid)
    
    # get whole sequence
    seq <- seqData(record)
        
    # get info on coi features
    coi <- coiInfo(record)
    
    # remove < and > from COI location span
    coi$loc <- gsub("<|>", "", coi$loc)
    
    # filter coi location info based on content
    # remove coi info if there is nesting of location operators
    coi$loc <- coi$loc[!grepl("\\(.*\\(.*", coi$loc)]
    # remove coi info if location operator order() is used
    coi$loc <- coi$loc[!grepl("order\\(", coi$loc)]
    
    # only proceed if there are remaining location spans
    if(length(coi$loc) > 0) {
      
      # condition 1:
      # coi location span follows format of digit..digit
      cond <- grepl("^[[:digit:]]*\\.\\.[[:digit:]]*$", coi$loc)
      if(any(cond)) {
        # restrict location spans 
        span <- coi$loc[cond]
        # take sequence start and end
        start <- as.numeric(gsub("\\.\\.[[:digit:]]*$", "", span))
        end <- as.numeric(gsub("^[[:digit:]]*\\.\\.", "", span))
        # cut sequence
        seqtemp <- seqCut(seq, start, end)
        # fill in f depending on length of seq
        for (k in 1:length(seqtemp)) { f <- c(f, id, seqtemp[k]) }
      }
      
      # condition 2:
      # location operator is "complement"
      cond <- grepl("^complement\\(", coi$loc)
      if(any(cond)) {
        # restrict location spans 
        span <- coi$loc[cond]
        # take sequence start and end
        start <- as.numeric(gsub("^complement\\(|\\.\\.[[:digit:]]*\\)$", "", span))
        end <- as.numeric(gsub("^complement\\([[:digit:]]*\\.\\.|\\)", "", span))
        # cut sequence
        seqtemp <- seqCut(seq, start, end)
        # check if there are unexpected character in a sequence
        x <- gsub("A|B|C|D|G|H|K|M|N|R|S|T|V|W|Y", "", seq)
        # print unexpected character if it is present
        if(nchar(x) > 0) { print(paste(x, "in", data$version)) } 
        # complement the sequence
        for(l in 1:length(seqtemp)) { 
          seqtemp[l] <- toupper(c2s(rev(comp(s2c(seqtemp[l]), ambiguous = TRUE)))) }
        # fill in f depending on length of seq
        for (k in 1:length(seqtemp)) { f <- c(f, id, seqtemp[k]) }  
      }
      
      # condition 3:
      # location operator is "join"
      cond <- grepl("^join\\(", coi$loc)
      if(any(cond)) {
        # restrict location spans 
        span <- coi$loc[cond]
        # take text within join operator
        span <- gsub("^join\\(|\\)$", "", span)
        # split span
        span <- strsplit(span,",")
        # join sequences per location span given
        for(m in 1:length(span)) {
          # take sequence start and end
          start <- as.numeric(gsub("\\.\\.[[:digit:]]*$", "", span[[m]]))
          end <- as.numeric(gsub("^[[:digit:]]*\\.\\.", "", span[[m]]))
          # cut and paste sequence
          seqtemp <- paste(seqCut(seq, start, end), collapse = "")
          f <- c(f, id, seqtemp)
        }
      } # end of condition 3
      
    } # if there are workable spans in the record
  } # per record in .gb file
  
  # write fasta if vector not empty
  if (length(f) > 0) { 
    # take file name and change to fasta
    temp <- gsub("_ncbi.gb", "_ncbi.fasta", files[i])
    
    # paste full path to file name (NCBI format)
    # temp <- paste("D:/Documents/NGS/entrez/wormstaxlist_gb/wormstaxlist_fasta_cut", "/", temp, sep = "")
    # paste full path to file name (kraken format)
    temp <- paste("D:/Documents/NGS/entrez/wormstaxlist_gb/wormstaxlist_fasta_cut_kraken", "/", temp, sep = "")
    
    write(f, temp) }
} # per .gb file

# inspect .fastas created
cut <- list.files(path = "D:/Documents/NGS/entrez/wormstaxlist_gb/wormstaxlist_fasta_cut", pattern = "_ncbi.fasta")
uncut <- list.files(path = "D:/Documents/NGS/entrez/wormstaxlist_gb/wormstaxlist_fasta", pattern = "_ncbi.fasta")

# list of sequences in uncut but not in cut
setdiff(uncut, cut)
locRaw[locRaw$file == "Polyplacotoma_mediterranea_ncbi.gb",]





# mine .gb files for taxid and organism 
# this is for use in validating taxids to be attached to BOLD sequences
# they will also be used to fill-in species with no available taxids from taxize (based on taxon or genbank accession)

# modify recInfo function to extract the following from each record:
# organism in addition to definition, version, taxid, source (range of sequence)
recInfo2 <- function(gbrecord) {
  # create empty list where to store record information
  list <- vector("list", length = 5)
  names(list) <- c("def", "version", "taxid", "source", "organism")
  
  # get accession version
  list$version <- gbrecord[grep("^VERSION", gbrecord)]
  list$version <- trimws(gsub("^VERSION", "", list$version))
  
  # get definition
  # line index containing start of definition
  ind1 <- grep("^DEFINITION", gbrecord)
  # indices of headers
  headers <- grep("^[[:alpha:]]", gbrecord)
  # line index of end of definition
  ind2 <- headers[which(headers == ind1) + 1] - 1
  list$def <- gbrecord[ind1:ind2]
  list$def <- paste(trimws(gsub("^DEFINITION", "", list$def)), collapse = " ")
  # remove last period
  list$def <- sub("[.]$", "", list$def)
  
  # line index containing "source"
  ind1 <- grep("^     source", gbrecord)
  # indices of feature starts
  headers <- grep("^     [[:alpha:]]", gbrecord)
  # line index of end of source key
  ind2 <- headers[which(headers %in% ind1) + 1] - 1
  
  # is there only one source key in the record?
  if (length(ind1) == 1) {
    # restrict actions to section under source key
    sub <- gbrecord[ind1:ind2]
    # get taxid
    list$taxid <- gsub("\\\"|\\s|/", "", sub[grep("db_xref=\"taxon:", sub)])
    list$taxid <- sub("db_xref=taxon:", "", list$taxid)
    # if there is no taxid, replace taxid with NA, and print message
    if (length(list$taxid) == 0) { 
      list$taxid <- NA
      print(paste("recInfo: sequence tagged", list$version, "has no taxid"))}
    
    # get location range of sequence
    list$source <- gsub("source|\\s", "", gbrecord[ind1])
    
    # get index of "organism" within sub
    ind <- grep("/organism=", sub)
    # if there is no organism, print message
    # if there is at least one index, take 1st index
    if(length(ind) == 0) { print(paste("recInfo2: there is no organism in sequence tagged", list$version)) } else { 
      list$organism <- trimws(gsub("\\\"|/|organism=", "", sub[ind[1]])) 
      # if there is more than one organism, print message
      if(length(ind) > 1) { print(paste("recInfo2: there are", length(ind), "organisms in sequence tagged", list$version))}
        }
  } else { print(paste("recInfo: there is more than one source key in sequence tagged", list$version))}
  
  # output list of record information 
  return(list)
}

# files to process
files <- list.files(pattern = "_ncbi.gb")
# non-fasta files in folder
list.files()[grep("_ncbi.gb", invert = TRUE, list.files())]

# run mining (1st run 2058.78 elapsed minutes)
for(i in 1:length(files)) {
  gblist <- gbParse(files[i])

  # do actions per record in the .gb file
  for (j in 1:length(gblist)) {
    # take jth record
    record <- gblist[[j]]
    # get sequence info
    data <- recInfo2(record)
    # add file name
    data$file <- files[i]

    # proceed only if there is one value per variable in list
    if(all(lengths(data) == rep(1, 6))) {
      data <- as.data.frame(data)
      write.table(data, file = "taxid_table_genbank.txt",
                  append = TRUE, row.names = FALSE, col.names = FALSE,
                  quote = TRUE, sep = "|", )
    } else { print(paste("Check sequence record tagged", list$version))}
  }
}

# read resulting table to add variable names
taxid_table_gb <- read.table(file = "taxid_table_genbank.txt", sep = "|", comment.char = "")
names(taxid_table_gb) <- c("def", "version",  "taxid",  "source", "organism", "file" )
# write table again
write.table(taxid_table_gb, file = "taxid_table_genbank.txt", row.names = FALSE, col.names = TRUE, quote = TRUE, sep = "|")





############################
# CONVERSION OF 16S GENBANK FILES TO FASTAS

# working directory: path of .gb files
setwd("D:/Documents/NGS/entrez/wormstaxlist_gb_16S/raw/")
# check number of sequences expected
# count after > is number after re-downloading bad files
# grep "^LOCUS" *_ncbi.gb -w | wc -l = 62794
# grep "^DEFINITION" *_ncbi.gb -w | wc -l = 62794
# grep "^ORIGIN" *_ncbi.gb -w | wc -l = 62727
# grep "^ORIGIN\|^WGS \|^CONTIG" *_ncbi.gb -w | wc -l = 62794
# grep "^//" *_ncbi.gb -w | wc -l = 62794
# grep "^     source" *_ncbi.gb -w | wc -l = 62796
# grep "db_xref=\"taxon:" -n *_ncbi.gb | wc -l = 62796

# list of .gbs to process = 4930
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
  m6 <- grep("^     source", gb) # there may be gbs with more that one "source"
  # compile in list
  occ <- list(m1, m2, m3, m4, m6)
  lens <- sapply(occ, length)
  
  # return file name and counts per marker if length of occurrences are not equal
  if (min(lens) == max(lens)) { NULL } 
  else { return(c(gbfile, lens, length(m5))) } 
}

# apply function over entire list of files
badfiles <- unlist(lapply(files, check))

# split vector
badmat <- as.data.frame(matrix(badfiles, byrow = TRUE, ncol = 7))
# format dataframe
colnames(badmat) <- c("file", "locus", "def", "ori", "end", "source", "wgs")
badmat[,-1] <- lapply(badmat[,-1], as.numeric)
# check counts
all(badmat$locus == badmat$def)
all(badmat$locus == badmat$end)
all(badmat$locus == (badmat$ori + badmat$wgs))
colSums(badmat[,-1])
# consistent with grep counts in Ubuntu:
# ORIGIN is not a consistent marker in each record
# WGS or CONTIG replace ORIGIN, but there are no sequences
badmat[badmat$locus != badmat$source,]
# Hordeum_vulgare_ncbi.gb and Pecten_maximus_ncbi.gb have gb records with more than 1 source
#                           file locus  def  ori  end source wgs
# 1            Bacillus__ncbi.gb  6128 6128 6108 6128   6128  20
# 2       Clunio_marinus_ncbi.gb     2    2    1    2      2   1
# 3     Corynebacterium__ncbi.gb   806  806  802  806    806   4
# 4     Escherichia_coli_ncbi.gb   250  250  241  250    250   9
# 5      Hordeum_vulgare_ncbi.gb     1    1    1    1      2   0
# 6       Mycobacterium__ncbi.gb   718  718  712  718    718   6
# 7       Pecten_maximus_ncbi.gb     7    7    7    7      8   0
# 8         Penicillium__ncbi.gb     3    3    2    3      3   1
# 9      Photobacterium__ncbi.gb   205  205  200  205    205   5
# 10  Propionibacterium__ncbi.gb   115  115   95  115    115  20
# 11 Salmonella_enterica_ncbi.gb    29   29   28   29     29   1

# some downloaded sequences do not contain 16S exclusively
# need to cut out 16S sequences only
# load following functions written previously for COI:
# 1) gbParse: takes .gb file and separates records
# 2) recInfo: takes record, returns list of definition, version, taxid, source (range of sequence)

# sxsInfo modified from coiInfo funciton
# 3) sxsInfo
# function to extract gene and location information of 16S
sxsInfo <- function(gbrecord) {
  # create empty vectors where locations and 16S gene names will be compiled
  loc <- vector()
  gene <- vector()
  
  # get accession version
  ver <- gbrecord[grep("^VERSION", gbrecord)]
  ver <- trimws(gsub("^VERSION", "", ver))
  
  # get features with "rRNA or gene"
  ind1 <- grep("^     rRNA|^     gene", gbrecord)
  # line indices where all features start + index of ORIGIN 
  headers <- grep("^     [[:alpha:]]|^ORIGIN", gbrecord)
  # indices where each "rRNA" feature ends
  ind2 <- headers[which(headers %in% ind1) + 1] - 1
  
  # only proceed if there is "^     rRNA or gene" in record
  if (length(ind1) > 0) {
    for (j in 1:length(ind1)) {
      # examines lines per gene feature
      sub <- gbrecord[ind1[j]:ind2[j]]
      # 16S pattern to search for 
      pattern <- "/gene.*16S"
      
      # if any 16S pattern is in the gene feature, extract location and gene name
      if (any(grepl(pattern, sub, ignore.case = TRUE))) { 
        # which lines contain the start of location
        stLine <- grep("^     rRNA|^     gene", sub)
        # which lines contain the end of location
        endLine <- grep("^                     /", sub)[1] - 1
        # location line
        l <- sub[stLine:endLine]
        # remove large spaces and collapse as one line
        l <- gsub("^                     ", "", l)
        l <- paste(l, collapse = "")
        
        # get location
        loc <- c(loc, gsub("gene|\\s|rRNA", "", l))
        # get gene
        gene <- c(gene, trimws(gsub("/gene=|\"", "", sub[grep("/gene=", sub)])))
      }
    }
    # if no COI pattern found in sections, print message in console
    if (length(gene) == 0) { 
      # print message
      print(paste("sxsInfo: no 16S gene information in sequence tagged", ver)) }
    
    # put gene and loc in dataframe to check duplicates
    d <- unique(data.frame(loc, gene))
    # re-assign back to vectors
    loc <- d$loc
    gene <- d$gene
    # return named list
    return(list(gene = gene, loc = loc))
  } else {
    # print message
    print(paste("sxsInfo: no gene information in sequence tagged", ver))
  }
}

# location data is not simply encoded, so need to check possible characters that arise in location info 
# all .gb files
files <- list.files(pattern = "_ncbi.gb")
# run gbParse, recInfo, sxsInfo (6.468833 mins elapsed time)
for (i in 1:length(files)) {
  # parse .gb files
  gblist <- gbParse(files[i])
  # nested loop
  for (j in 1:length(gblist)) {
    # perform actions per record
    record <- gblist[[j]]
    # get record and gene info
    l1 <- recInfo(record)
    l2 <- sxsInfo(record)
    
    # base nrows on number of genes listed as 16S
    df <- data.frame(gene = l2$gene, loc = l2$loc)
    # add unique data per record
    df$def <- l1$def
    df$version <- l1$version
    df$taxid <- ifelse(is.null(l1$taxid), NA, l1$taxid) 
    df$source <- ifelse(is.null(l1$source), NA, l1$source)
    df$file <- files[i]
    
    # append df to file
    write.table(df, "D:/Documents/NGS/entrez/wormstaxlist_gb_16S/locInfo_16S.txt", append = TRUE, col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  }
}
# console outputs consistent with badmat

# check location information 
locRaw <- read.table("D:/Documents/NGS/entrez/wormstaxlist_gb_16S/locInfo_16S.txt", sep = "\t", quote = "", comment.char = "") # 63050 rows
# add names to columns
names(locRaw) <- c("gene", "loc", "def", "version", "taxid", "source", "file")

# check the format of source
# confirm if it always starts with 1..(length of sequence)
temp <- locRaw$source
# do all begin with 1..
all(grepl("^1\\.\\.", temp))
all(grepl("^1\\.\\.", temp[!is.na(temp)]))
locRaw[!grepl("^1\\.\\.", temp), ]
# do all follow the format: 1..digit
all(grepl("^1\\.\\.[[:digit:]]*$", temp[!is.na(temp)]))
unique(sub("^1\\.\\.[[:digit:]]*$", "", temp[!is.na(temp)]))
# sequence lengths
len <- as.integer(sub("^1\\.\\.", "", temp[!is.na(temp)]))
# range of sequence lengths: 46 to 5136383
range(len)
# distribution of lengths
hist(len)
hist(sort(len)[1:61000], breaks = 100)
# 96.90% of lengths are less than 1600 bp
sum(len < 1600) / length(len)

# which accession #s are repeated?
df <- data.frame(table(locRaw$version))
# vector of accession numbers that appear more than once
dups <- df[df$Freq > 1,]$Var1
# print duplicates based on accession #
temp <- locRaw[locRaw$version %in% dups,]
temp <- temp[order(temp$version),]
head(temp[, c("version", "taxid", "source", "file")])
# some duplicates: Bacillus__ncbi.gb vs. Bacillus_sp._ncbi.gb

# number of unique accession #s in df showing duplicates
length(unique(temp$version)) # 11252
# number of unique entries based on record contents alone
nrow(unique(temp[,c('def', 'version', 'taxid', 'source')]))
# note: number of version is indicative of the #s of unique records

# remove duplicates due to species synonyms to see other duplicates 
# duplication will then be due to multiple 16S locations per record
# remove file column
loc <- locRaw[, -7]
# remove duplicates based on all remaining variables
loc <- loc[!duplicated(loc), ]
# check again which accession numbers are repeated
df <- data.frame(table(loc$version))
# vector of accession numbers that appear more than once
dups <- as.character(df[df$Freq > 1,]$Var1)
# print duplicates based on accession #
temp <- loc[loc$version %in% dups,]
temp <- temp[order(temp$version),]
temp[, c("version", "source", "loc", "gene")]
# there are 186 records associated with more than one 16S gene location 
length(unique(temp$version))

# sxsChecker
# function for printing all 16S features in .gb record if there is more than one
sxsChecker <- function(gbrecord) {
  # create empty vectors where locations and 16S gene names will be compiled
  loc <- vector()
  gene <- vector()
  # vector where feature section will be stored
  feature <- vector()
  
  # get features with "gene" or "rRNA"
  ind1 <- grep("^     gene|^     rRNA", gbrecord)
  # line indices where all features start + index of ORIGIN 
  headers <- grep("^     [[:alpha:]]|^ORIGIN", gbrecord)
  # indices where each "gene"or "rRNA" feature ends
  ind2 <- headers[which(headers %in% ind1) + 1] - 1
  
  for (j in 1:length(ind1)) {
    # examines lines per gene feature
    sub <- gbrecord[ind1[j]:ind2[j]]
    # 16S pattern to search for 
    pattern <- "/gene.*16S"
    
    # if any 16S pattern is in the gene feature, extract location and gene name
    if (any(grepl(pattern, sub, ignore.case = TRUE))) { 
      loc <- c(loc, gsub("gene|\\s|rRNA", "", sub[grep("^     gene|^     rRNA", sub)]))
      gene <- c(gene, trimws(gsub("/gene=|\"", "", sub[grep("/gene=", sub)])))
      # save feature in vector
      feature <- c(feature, sub)
    }
  }
  
  # get accession version
  ver <- gbrecord[grep("^VERSION", gbrecord)]
  ver <- trimws(gsub("^VERSION", "", ver))
  
  # if no COI pattern found in sections, print message in console
  if (length(gene) == 0) { 
    # print message
    print(paste("geneInfo: no 16S gene information in sequence tagged", ver)) }
  
  # if there is more than one 16S feature, print accession # and feature
  if (length(gene) > 1) {
    print(ver)
    print(feature)
  }
}

# vector of accessions associated with more than one 16S location
vers <- unique(temp$version)
# find file names linked to example accession #s
f <- unique(locRaw[locRaw$version %in% vers, ]$file) 
# mapping of accession version with file name
unique(locRaw[locRaw$version %in% vers, c("version", "file")])

# print sections when 16S locations > 1
con <- file("D:/Documents/NGS/entrez/wormstaxlist_gb_16S/sxsChecker_log.txt")
sink(con, append=TRUE)
# loop
for (i in 1:length(f)) {
  gblist <- gbParse(f[i])
  
  for (j in 1:length(gblist)) {
    # perform actions per record
    record <- gblist[[j]]
    sxsChecker(record)
  }
}
sink()

# inspect sequence spans with location operators
d <- unique(locRaw$loc)
d[grep("[[:alpha:]]", d)] # 115 / 3792 or 3%
# no sequence spans with  nested location operators
d[grep("\\(.*\\(", d)]

# inspect sequence spans with no location operators
d <- unique(locRaw$loc)
d <- d[grep("[[:alpha:]]", d, invert = TRUE)]
# removing < and >, do all sequence spans follow the format of digit..digit
t <- gsub("<|>", "", d)
all(grepl("^[[:digit:]]*\\.\\.[[:digit:]]*$", t))
# range of sequence lengths
start <- as.integer(sub("\\.\\.[[:digit:]].*", "", t))
end <- as.integer(sub("[[:digit:]].*\\.\\.", "", t))
len <- end - start
hist(len, breaks = 50)
range(len) 
# find sequences with length less than 100 bp
smallspans <- d[which(len < 100)]
locRaw[locRaw$loc %in% smallspans,-c(3,5)]
# some sequences are short since the record is an amalgamation of different regions.
# e.g. HE612855.1 Acanthaster planci mitochondrial partial D-loop and 16S rRNA gene, isolate KIN21. whose 16S span is 460..>509
# there are also sequences that are less than 10 bp:
locRaw[locRaw$loc %in% d[which(len < 10)],-c(3,5)]
# e.g.: AJ243775.1 Cordyceps sinensis partial 16S rRNA gene, 5.8S rRNA gene, partial 26S rRNA; ITS1 and 2, sclerotium. with 16S span <1..4, but source span is 1..588

# % of records whose 16S location span is less than 100 = 1.513% = 954/63050
(nrow(locRaw[locRaw$loc %in% smallspans,]) / nrow(locRaw))*100
# % of sequences with small spans whose location span == source span; 17/954
t <- locRaw[locRaw$loc %in% smallspans,]
t$loc <- gsub("<|>", "", t$loc)
t[which(t$loc == t$source), -c(3,5)]
# % of sequences whose COI span == source span is 89.45%
t <- locRaw
t$loc <- gsub("<|>", "", t$loc)
(nrow(t[t$loc == t$source,]) / nrow(t))*100

# load needed functions:
# 1) seqData: function for extracting sequence from GenBank record
# 2) seqCut: function for cutting a given sequence based on a location span

# files for conversion to 16S fastas
files <- list.files(pattern = "_ncbi.gb")

# file which produced error because of single base position
# record <- gbParse("Bacillus__ncbi.gb")
# record <- record[[5981]]

# for loop for converting .gbs to .fastas
# load required package:
library(seqinr)
# sequences are restricted to COI locations
# 7.717333 elapsed minutes to run (kraken format) 
for (i in 1:length(files)) {
  gblist <- gbParse(files[i])
  
  # empty vector where fasta will be assembled per .gb
  f <- vector()
  # do actions per record in the .gb file
  for (j in 1:length(gblist)) {
    # take jth record
    record <- gblist[[j]]
    
    # get record metadata 
    data <- recInfo(record)
    
    # id info for sequence (unique for each record)
    # id <- paste0(">", data$version, " ", data$def)
    
    # id if NCBI format will be overriden for kraken format
    id <- paste0(">", "kraken:taxid|", data$taxid)
    
    # get whole sequence
    seq <- seqData(record)
    
    # get info on 16S features
    sx <- sxsInfo(record)
    
    # remove < and > from 16S location span
    sx$loc <- gsub("<|>", "", sx$loc)
    
    # filter 16S location info based on content
    # remove 16S info if there is nesting of location operators
    sx$loc <- sx$loc[!grepl("\\(.*\\(.*", sx$loc)]
    # remove 16S info if location operator order() is used
    sx$loc <- sx$loc[!grepl("order\\(", sx$loc)]
    
    # only proceed if there are remaining location spans
    if(length(sx$loc) > 0) {
      
      # condition 1:
      # 16S location span follows format of digit..digit
      cond <- grepl("^[[:digit:]]*\\.\\.[[:digit:]]*$", sx$loc)
      if(any(cond)) {
        # restrict location spans 
        span <- sx$loc[cond]
        # take sequence start and end
        start <- as.numeric(gsub("\\.\\.[[:digit:]]*$", "", span))
        end <- as.numeric(gsub("^[[:digit:]]*\\.\\.", "", span))
        # cut sequence
        seqtemp <- seqCut(seq, start, end)
        # fill in f depending on length of seq
        for (k in 1:length(seqtemp)) { f <- c(f, id, seqtemp[k]) }
      }
      
      # condition 2:
      # location operator is "complement"
      cond <- grepl("^complement\\(", sx$loc) & grepl("\\.\\.", sx$loc)
      if(any(cond)) {
        # restrict location spans 
        span <- sx$loc[cond] 
        # take sequence start and end and catch errors
        # tryCatch({start <- as.numeric(gsub("^complement\\(|\\.\\.[[:digit:]]*\\)$", "", span))}, warning=function(w) print(c(data$version, j)))
        # tryCatch({end <- as.numeric(gsub("^complement\\([[:digit:]]*\\.\\.|\\)", "", span))}, warning=function(w) print(c(data$version, j)))
        start <- as.numeric(gsub("^complement\\(|\\.\\.[[:digit:]]*\\)$", "", span))
        end <- as.numeric(gsub("^complement\\([[:digit:]]*\\.\\.|\\)", "", span))
        # cut sequence
        seqtemp <- seqCut(seq, start, end)
        # check if there are unexpected character in a sequence 
        x <- gsub("A|B|C|D|G|H|K|M|N|R|S|T|V|W|Y", "", seq)
        # print unexpected character if it is present
        if(nchar(x) > 0) { print(paste(x, "in", data$version)) } 
        # complement the sequence
        for(l in 1:length(seqtemp)) { 
          seqtemp[l] <- toupper(c2s(rev(comp(s2c(seqtemp[l]), ambiguous = TRUE)))) }
        # fill in f depending on length of seq
        for (k in 1:length(seqtemp)) { f <- c(f, id, seqtemp[k]) }  
      }
      
      # condition 3:
      # location operator is "join"
      cond <- grepl("^join\\(", sx$loc)
      if(any(cond)) {
        # restrict location spans 
        span <- sx$loc[cond]
        # take text within join operator
        span <- gsub("^join\\(|\\)$", "", span)
        # split span
        span <- strsplit(span,",")
        # join sequences per location span given 
        for(m in 1:length(span)) {
          # take sequence start and end
          start <- as.numeric(gsub("\\.\\.[[:digit:]]*$", "", span[[m]]))
          end <- as.numeric(gsub("^[[:digit:]]*\\.\\.", "", span[[m]]))
          # cut and paste sequence
          seqtemp <- paste(seqCut(seq, start, end), collapse = "")
          f <- c(f, id, seqtemp)
        }
      } # end of condition 3
      
    } # if there are workable spans in the record
  } # per record in .gb file
  
  # write fasta if vector not empty
  if (length(f) > 0) { 
    # take file name and change to fasta
    temp <- gsub("_ncbi.gb", "_ncbi.fasta", files[i])
    
    # paste full path to file name (NCBI format)
    # temp <- paste("D:/Documents/NGS/entrez/wormstaxlist_gb_16S/wormstaxlist_fasta_cut", "/", temp, sep = "")
    # paste full path to file name (kraken format)
    temp <- paste("D:/Documents/NGS/entrez/wormstaxlist_gb_16S/fasta_cut_kraken", "/", temp, sep = "")
    
    write(f, temp) }
} # per .gb file





##### run codes here:
ptm <- proc.time()
# put function here:

proc.time() - ptm
beepr::beep(5)