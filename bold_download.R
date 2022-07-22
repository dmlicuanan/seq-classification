# 17 July 2022
# download sequences from BOLD
# https://docs.ropensci.org/bold/#large-data

# install packages
install.packages("bold")
install.packages("taxize")
# load packages
library(bold)
library(taxize)
library(stringr)

# taxize authenthication -- run everytime
options(ENTREZ_KEY = "7d21c9bc180e2c7f4c118b010e6f236f3008")

# set working directory
setwd("D:/Documents/NGS/entrez/")
# read-in taxa list
wormstaxlist <- readLines("wormstaxlist")
# arrange alphabetically
wormstaxlist <- wormstaxlist[order(wormstaxlist)]
# remove "'" since it causes Internal Server Error (HTTP 500) in BOLD
wormstaxlist <- gsub("\\'", "", wormstaxlist)
# check other symbols/punctuations in list
wormstaxlist[grep("[[:punct:]]", wormstaxlist)]

# create directory where fastas will be saved
# dir.create("./bold_fastas")
setwd("./bold_fastas")

# function for getting fasta (modified from Mares)
get_fasta <- function(taxon, filename) {
  x <- bold_seqspec(taxon = taxon, marker = c("COI-5P", "COI-3P"))
  x <- x[x$markercode == "COI-5P" | x$markercode == "COI-3P", ]
  x[x==""] <- NA 
  b_acc <- x$processid
  b_tax <- ifelse(!is.na(x$species_name), x$species_name, ifelse(!is.na(x$genus_name),x$genus_name, ifelse(
    !is.na(x$family_name), x$family_name,ifelse (
      !is.na(x$order_name), x$order_name,ifelse(
        !is.na(x$class_name), x$class_name,x$phylum_name)))))
  b_mark <- x$markercode
  n_acc <- ifelse(!is.na(x$genbank_accession),ifelse(!is.na(x$genbank_accession),paste0("|",x$genbank_accession),""),"")
  
  seq <- x$nucleotides
  seqname <- paste(b_acc, b_tax, b_mark, sep="|")
  seqname<-paste0(seqname, n_acc)
  Y <- cbind(seqname, seq)
  colnames(Y)<-c("name","seq")
  fastaLines = c()
  for (rowNum in 1:nrow(Y)){
    fastaLines = c(fastaLines, as.character(paste(">", Y[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines, as.character(Y[rowNum,"seq"]))
  }
  writeLines(fastaLines, filename)
}

# assign list of taxa
list <- wormstaxlist

sink("log.txt", append = TRUE)
# for-loop length(list)
for(i in 24991:length(list)) {
  # check if taxon in BOLD before proceeding
  taxon <- list[i]
  x <- bold_tax_name(taxon)
  
  # only proceed if taxon is in BOLD 
  if(dim(x)[2] > 2) {
    # progress
    cat("Processing", taxon, ":", i, "of", length(list), "\n")
    
    tryCatch({
      get_fasta(taxon = taxon, filename = paste0(gsub(" ", "_", taxon), "_bold.fasta"))
    }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  } else { NULL }
}
sink()

# Error in curl::curl_fetch_memory(x$url$url, handle = x$url$handle) : 
# Timeout was reached: [v4.boldsystems.org] Operation timed out after # 10015 milliseconds with 0 out of 0 bytes 

# testing if can record changes to files.

# will remove "received in line 80
