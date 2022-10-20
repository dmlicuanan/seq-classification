# creates a list of species using a WoRMS dataset

# load packages
library(readr)
library(stringr)

# use file requested from WoRMS site: https://www.marinespecies.org/download/
# WoRMS downloads: This page provides access to a monthly copy of WoRMS in DwC-A format. Only taxonomic data is provided, extra data fields may be available upon request. Download access is limited to 1 year, starting from the last download. Data originating from AlgaeBase is excluded, since the license does not allow redistribution. Old (legacy) copies before 2017, in MS Access format, are still available here.

# set working directory
setwd("D:/Documents/NGS/input/WoRMS_download_2022-07-01")

# preview file
readLines("taxon.txt", 10)
readLines("speciesprofile.txt", 10)
readLines("vernacularname.txt", 10)

# read in file
raw <- readr::read_delim("taxon.txt", delim = "\t")
# remove some NAs
raw <- as.data.frame(raw[which(!is.na(raw$acceptedNameUsage)),])
head(raw)

# explore variables
colnames(raw)
table(raw$taxonRank)
table(raw$taxonomicStatus)

# filtering
# 1) restrict to rows at Species rank, and whose taxonomic status is accepted
d <- raw[which(raw$taxonRank == "Species" & raw$taxonomicStatus == "accepted"), ]
# 2) remove prokaryotes
proks <- c("Bacteria", "Archaea", "Viruses")
d <- d[which(!(d$kingdom %in% proks)),]
# 3) remove rows with NA in genus and specificEpithet columns
d <- d[which(!(is.na(d$genus))|!(is.na(d$specificEpithet))),]

# check rows where scientificName is not the same as acceptedNameUsage
d[which(d$acceptedNameUsage != d$scientificName),]
# just one entry: everything else have scientificName the same as accpetedNameUsage

# check word counts for acceptedNameUsage or scientificName
table(str_count(d$scientificName, "\\w+"))
d[which(str_count(d$scientificName, "\\w+") >= 3), c(6,7,16,18)]

# to simplify search term, we just paste genus with specificEpithet
# https://dwc.tdwg.org/terms
d$term <- paste(d$genus, d$specificEpithet)
length(d$term)
length(unique(d$term))
# unique values only:
wormstaxlist <- d$term

# which have special characters in string
wormstaxlist[grep("[[:punct:]]", wormstaxlist)]

# write file readable in Unix
f <- file("D:/Documents/NGS/entrez/wormstaxlist", open="wb")
cat(wormstaxlist, file = f, sep = "\n")
close(f)





# wormstaxlist is incomplete since there are 52621 taxids not represented in the worms_coi.fasta
# for COI sequences, we just extracted the corresponding sequences in the MARES database then added them to worms_coi.fasta (see bold_format.R)
# for 16S sequences, the sequences need to be downloaded from NCBI, so we need a taxa list 

# some codes below from bold_format.R
# readFasta function
readFasta <- function(file){
  raw <- readLines(file)  
  indNames <- grep(">", raw)  
  stSeq <- indNames + 1  
  enSeq <- c(indNames[-1] - 1, length(raw))  
  seqOut <- mapply(function(s, e, inSeq){paste(inSeq[s:e], collapse = "")},
                   s = stSeq, e = enSeq, MoreArgs = list(inSeq = raw))  
  id <- gsub(">", "", raw[indNames])  
  return(data.frame(id = id, seq = seqOut, stringsAsFactors = FALSE))
}

# read in MARES database fasta
mares <- readFasta("D:/Documents/NGS/mares/mares_nobar_taxonomy/MARES_NOBAR_BOLD_NCBI_sl_kraken.fasta")
# add taxid column
mares$taxid <- gsub("kraken:taxid\\|", "", mares$id)

# read in worms fasta database
worms <- readFasta("D:/Documents/NGS/entrez/kraken_fastas/worms_coi_deduplicated.fasta")
# add taxid column
worms$taxid <- gsub("kraken:taxid\\|", "", worms$id)

# taxids in MARES but not in worms fasta
miss <- setdiff(unique(mares$taxid), unique(worms$taxid))

# read in nodes.dmp from NCBI taxdump to get corresponding taxa name of taxids
names <- data.table::fread("D:/Documents/NGS/entrez/taxonomy/names.dmp", sep = "|", drop = 5, col.names = c("taxid", "name", "name_unique", "name_class"))
# variables are:
# tax_id = the id of node associated with this name
# name_txt = name itself
# unique name = the unique variant of this name if name not unique
# name class = (synonym, common name, ...)

# remove tabs from all columns
names <- as.data.frame(lapply(names, trimws))

# reduce names so that only scientific names appear
names <- names[names$name_class == "scientific name", ]

# compile results
df <- data.frame(taxid = miss, 
                 taxon = names$name[match(miss, names$taxid)],
                 name_unique = names$name_unique[match(miss, names$taxid)])

# load library
library(stringr)

# get taxa that are not empty/NA
taxa <- df$taxon[!is.na(df$taxon) & df$taxon != ""]
# remove parentheses and all text within
taxa <- str_replace(taxa, " \\s*\\([^\\)]+\\)", "")
# remove NAs and get unique values
taxa <- unique(taxa[!is.na(taxa)])

# clean 1st and 2nd words separately
temp <- data.frame(genus = word(taxa, 1), species = word(taxa, 2))

# for 1st word:
# only keep taxa if genus starts with capital letter
temp <- temp[grep("^[A-Z]", temp$genus),]
# remove taxa if there is ' - ( ) in genus
temp <- temp[grep("[\\'\\-\\(\\)]", temp$genus, perl = TRUE, invert = TRUE), ]

# for 2nd word:
# replace second word with "" if starts with a capital letter
temp$species <- ifelse(grepl("^[A-Z]", temp$species), "", temp$species)
# replace NA with ""
temp$species <- ifelse(is.na(temp$species), "", temp$species)
# remove rows with quotes in species epithet
temp <- temp[grep("[\\']", temp$species, perl = TRUE, invert = TRUE), ]

# empty vector for collapsed genus + species
taxa <- vector(mode = "character", nrow(temp))
# collapse strings together
for(i in 1:nrow(temp)) { 
  taxa[i] <- trimws(paste(temp$genus[i], temp$species[i], collapse = " ")) }

# remove taxa with unexpected symbols:
taxa[grep(pattern = "^[a-zA-Z -.]*$", invert = TRUE, taxa)]
taxa <- unique(taxa[grep(pattern = "^[a-zA-Z -.]*$", taxa)])

# check which have special characters in string, besides period
s <- gsub("\\.", "", taxa)
s[grep("[[:punct:]]", s)]

# write file readable in Unix
f <- file("D:/Documents/NGS/entrez/marestaxlist_uniq", open="wb")
cat(taxa, file = f, sep = "\n")
close(f)