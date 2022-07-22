# this script is for downloading fasta files from BOLD
# gene: COI 
# species: marine eukaryotes
# see https://docs.ropensci.org/bold/#large-data 

# install packages
install.packages("bold")
# load packages
library(bold)

# taxize authenthication -- run everytime: but needed only for use of taxize package
# API key from amlicuanan account at https://www.ncbi.nlm.nih.gov/account/settings/
# options(ENTREZ_KEY = "7d21c9bc180e2c7f4c118b010e6f236f3008")

# set working directory to where wormstaxlist is located
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

# 1) for laptop run (1st half of species list)
sink("log.txt", append = TRUE)
for(i in 24991:111160) {
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
# Timeout was reached: [v4.boldsystems.org] Operation timed out after # 10015 milliseconds with 0 out of 0 bytes received

# 2) for 172 PC run (start at 2nd half of species list)
sink("log_172.txt", append = TRUE)
for(i in 111161:length(list)) {
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
