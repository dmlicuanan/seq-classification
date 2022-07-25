# this script is for downloading fasta files from BOLD
# gene: COI 
# species: marine eukaryotes; species list used contained in "wormstaxlist" file
# see https://docs.ropensci.org/bold/#large-data 

# install packages
# install.packages("bold")
# load packages
library(bold)

# taxize authenthication -- run everytime: but needed only for use of taxize package
# API key from amlicuanan account at https://www.ncbi.nlm.nih.gov/account/settings/
# options(ENTREZ_KEY = "7d21c9bc180e2c7f4c118b010e6f236f3008")

# set working directory to where wormstaxlist is located 
setwd("D:/Documents/NGS/entrez/")
# read-in taxa list
wormstaxlist <- readLines("wormstaxlist", encoding = "UTF-8")
# arrange alphabetically
wormstaxlist <- wormstaxlist[order(wormstaxlist)]
# remove "'" since it causes Internal Server Error (HTTP 500) in BOLD
wormstaxlist <- gsub("\\'", "", wormstaxlist)
# check other symbols/punctuations in list
# wormstaxlist[grep("[[:punct:]]", wormstaxlist)]

# assign list of taxa
list <- wormstaxlist
rm(wormstaxlist)

# create directory where fastas will be saved
# dir.create("./bold_fastas")
setwd("./bold_fastas")

# note: 
# a for-loop containing a custom function (get_fasta) was used to download BOLD sequences
# both for-loop and get_fasta were modified during the run to remove errors
# version 1 of codes (lines commented below) were used for species approx. indexed 1:27203 and 111161:112648
# 1:27203 was run in laptop; 111161:112648 run in 172 PC remotely
# version 2 of codes were used for the rest of the sequences
# laptop run finished early, tail of wormstaxlist (around 215000:222320) was run in laptop instead of 172 PC

#################################
# version 1 of get_fasta function
# modified from https://github.com/wpearman1996/MARES_database_pipeline/blob/master/step2_retrieve_bold.r
# get_fasta <- function(taxon, filename) {
#   x <- bold_seqspec(taxon = taxon, marker = c("COI-5P", "COI-3P"))
#   x <- x[x$markercode == "COI-5P" | x$markercode == "COI-3P", ]
#   x[x==""] <- NA 
#   b_acc <- x$processid
#   b_tax <- ifelse(!is.na(x$species_name), x$species_name, ifelse(!is.na(x$genus_name),x$genus_name, ifelse(
#     !is.na(x$family_name), x$family_name,ifelse (
#       !is.na(x$order_name), x$order_name,ifelse(
#         !is.na(x$class_name), x$class_name,x$phylum_name)))))
#   b_mark <- x$markercode
#   n_acc <- ifelse(!is.na(x$genbank_accession),ifelse(!is.na(x$genbank_accession),paste0("|",x$genbank_accession),""),"")
#   
#   seq <- x$nucleotides
#   seqname <- paste(b_acc, b_tax, b_mark, sep="|")
#   seqname<-paste0(seqname, n_acc)
#   Y <- cbind(seqname, seq)
#   colnames(Y)<-c("name","seq")
#   fastaLines = c()
#   for (rowNum in 1:nrow(Y)){
#     fastaLines = c(fastaLines, as.character(paste(">", Y[rowNum,"name"], sep = "")))
#     fastaLines = c(fastaLines, as.character(Y[rowNum,"seq"]))
#   }
#   writeLines(fastaLines, filename)
# }

# version 1 of for-loop
# log.txt for laptop run; log_172.txt for 172 PC run
# sink("log.txt", append = TRUE) 
# for(i in 1:length(list)) {
#   # check if taxon in BOLD before proceeding
#   taxon <- list[i]
#   x <- bold_tax_name(taxon)
#   
#   # only proceed if taxon is in BOLD 
#   if(dim(x)[2] > 2) {
#     # progress
#     cat("Processing", taxon, ":", i, "of", length(list), "\n")
#     
#     tryCatch({
#       get_fasta(taxon = taxon, filename = paste0(gsub(" ", "_", taxon), "_bold.fasta"))
#     }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
#   } else { NULL }
# }
# sink()
#################################

# version 2 of get_fasta function
# removed marker argument in bold_seqspec and stopped process if object is atomic (just NA) and marker is not COI-5P or COI-3P
# https://github.com/ropensci/bold/issues/76
get_fasta <- function(taxon, filename) {
  x <- bold_seqspec(taxon = taxon)
  
  if(!is.atomic(x) & ("COI-5P" %in% x$markercode | "COI-3P" %in% x$markercode)) {
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
  } else { NULL }
}

# version 2 of for-loop
# log.txt for laptop run; log_172.txt for 172 PC run
sink("log.txt", append = TRUE)
for(i in 1:length(list)) {
  # check if taxon in BOLD before proceeding
  taxon <- list[i]
  x <- bold_tax_name(taxon)
  s <- bold_stats(taxon)
  
  # only proceed if taxon is in BOLD and there are records
  if(dim(x)[2] > 2 & s$total_records > 1) {
    # progress
    cat("Processing", taxon, ":", i, "of", length(list), "\n")
    
    tryCatch({
      get_fasta(taxon = taxon, filename = paste0(gsub(" ", "_", taxon), "_bold.fasta"))
    }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  } else { NULL }
}
sink()

