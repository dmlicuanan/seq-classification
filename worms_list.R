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