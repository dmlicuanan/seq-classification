# this script is for downloading fasta files from BOLD
# gene: COI 
# species: marine eukaryotes; species list used contained in "wormstaxlist" file
# see https://docs.ropensci.org/bold/#large-data 

# install packages
# install.packages("bold")
# load packages
library(bold)

# taxize authenthication -- run everytime: but needed only for use of taxize package
# API key from NCBI account at https://www.ncbi.nlm.nih.gov/account/settings/
# options(ENTREZ_KEY = "API_KEY")

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
# log_rerun.txt for after double-checking errors
# sink("log_rerun.txt", append = TRUE)
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





# check errors in log.txt and log_172.txt 
# to see if some species fastas were missed due to errors
# bold package is needed, as well as version 2 of get_fasta above

# set working directory to path of log files
setwd("D:/Documents/NGS/entrez/bold_fastas")

# merge log.txt and log_172.txt in one object
log <- c(readLines("log.txt"), readLines("log_172.txt"))

# find all occurrences of errors
unique(log[grep("^Processing", log, invert = TRUE)])
# [1] "ERROR : BOLD servers returned an error - we're not sure what happened"                
# [2] " try a smaller query - or open an issue and we'll try to help "                       
# [3] "ERROR : $ operator is invalid for atomic vectors "                                    
# [4] "ERROR : Operation was aborted by an application callback "                            
# [5] "ERROR : Send failure: Connection was reset "                                          
# [6] "ERROR : Failure when receiving data from the peer "                                   
# [7] "ERROR : incomplete final line found by readTableHeader on 'text' "                    
# [8] "ERROR : subscript out of bounds "                                                     
# [9] "ERROR : Bad Gateway (HTTP 502) "                                                      
# [10] "ERROR : Timeout was reached: [v4.boldsystems.org] Send failure: Connection was reset "

# get species succeeded by "ERROR" messages
# indices of lines before error messages:
err <- grep("ERROR", log) - 1
# check if all contain a species name
all(grepl("Processing", log[err]))

# extract species for re-run/double-checking; contains 1205 species
list <- gsub("^Processing (.+) :.*", "\\1", log[err])

# re-run version 2 of get_fasta function above
# re-run version 2 of for-loop above, but use "log_rerun.txt" as sink
# index of species printed in "log_rerun.txt" is not the same as their index in wormstaxlist

# checking of log_rerun.txt
log <- readLines("log_rerun.txt")

# only two additional fastas were downloaded: 
# Peltogasterella_gracilis_bold.fasta and "Parapterois_heterura_bold.fasta"
# this means that not all species printed in log files have corresponding fastas
# for example: 
taxon <- gsub("^Processing (.+) :.*", "\\1", log[1])
# it satisfies if statement in for-loop -- hence it is printed in log file
x <- bold_tax_name(taxon)
s <- bold_stats(taxon)
dim(x)[2] > 2 & s$total_records > 1
# but it does not satisfy if statement in get_fasta, which comes after cat command:
x <- bold_seqspec(taxon = taxon)
!is.atomic(x) & ("COI-5P" %in% x$markercode | "COI-3P" %in% x$markercode)

# checking of log_rerun.txt for further errors
unique(log[grep("^Processing", log, invert = TRUE)])
# [1] "ERROR : incomplete final line found by readTableHeader on 'text' "

# get species succeeded by "ERROR" messages
err <- grep("ERROR", log) - 1
# extract species for which error occurred
taxa <- gsub("^Processing (.+) :.*", "\\1", log[err])
taxon <- taxa[2]

# find where error occurs: 
x <- bold_tax_name(taxon)
s <- bold_stats(taxon)
dim(x)[2] > 2 & s$total_records > 1
# error occurs here for both species, hence no need for further re-run:
x <- bold_seqspec(taxon = taxon)

# final notes: 
# format of BOLD fastas are as follows:
# >processid|species|markercode|genbank_accession
# processid = BOLD Process IDs are unique codes automatically generated for each new record added to a project. They serve to connect specimen information, such as taxonomy, collection data and images, to the DNA barcode sequence for that specimen.

# sample:
# >GBMIN123652-17|Parapterois heterura|COI-5P|KP267594
# CCTGTATCTGGTATTTGGTGCCTGAGCCGGCATAGTAGGCACAGCCTTAAGCCTACTTATTCGAGCAGAGCTTAGTCAACCAGGCGCTCTATTGGGAGACGACCAAATTTACAATGTAATTGTTACAGCACACGCTTTTGTAATAATTTTCTTTATAGTAATACCAATTATGATTGGGGGATTTGGAAATTGACTTATCCCACTAATGATCGGAGCCCCAGACATGGCATTCCCCCGAATAAATAATATGAGCTTCTGACTCTTGCCGCCCTCTTTCCTCCTCCTCCTTGCCTCTTCAGGTGTTGAAGCAGGAGCCGGAACAGGATGAACCGTTTACCCGCCCCTAGCGGGTAATCTTGCCCACGCAGGAGCATCCGTAGATTTAACAATCTTCTCTCTCCATTTAGCCGGAATTTCATCGATTCTGGGGGCAATCAACTTCATTACAACAATCATTAACATGAAACCCCCAGCGATTTCCCAGTACCAAACACCTCTATTTGTGTGGGCTGTACTGATTACTGCGGTACTTTTACTTCTCTCACTTCCAGTCCTCGCCGCTGGCATTACAATACTACTCACAGATCGAAACCTTAATACTACTTTCTTCGACCCGGCGGGAGGGGGAGACCCAATTCTGTACCAACACCTTTTC

# total number of fastas = 22350
length(list.files(pattern = "\\.fasta$"))
