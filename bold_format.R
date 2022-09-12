# script for formatting BOLD .fastas to kraken format
# ID string of BOLD .fastas is as follows:
# process ID|taxon|marker code|genbank accession (if present)
# examples:
# "GBSP4661-12|Acanthella acuta|COI-5P|HQ379408" 
# "SIMP208-19|Acanthastrea echinata|COI-5P"  

# custom functions:
# 1) readFasta()
# feads .fasta file and returns data frame with variables id and seq
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

# 2) writeFasta()
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

# set working directory to location of BOLD .fastas
setwd("D:/Documents/NGS/entrez/bold_fastas")
# non-fasta files in folder
list.files()[grep("_bold.fasta", invert = TRUE, list.files())]
# BOLD fastas to process:
files <- list.files(pattern = "_bold.fasta")

# extract ID lines and write to file
for(i in 1:length(files)) {
  # read fasta
  temp <- readFasta(files[i])
  # write ID to file
  write(temp$id, file = "bold_id_list", append = TRUE)
}

# read in file with BOLD ids
df <- read.table("bold_id_list", sep = "|", fill = TRUE)
# add names
names(df) <- c("processid", "taxon", "markercode", "genbank_accession")

# load required packages
library(taxize)

# species list that need taxids
speclist <- unique(df$taxon)
# generate table for taxids (elspased time: 20999.34)
for(i in 1:length(speclist)) {
  # get taxize::get_uid_ result
  res <- get_uid_(sci_com = speclist[i], 
                  key = "81af93ade140a45a355b621d22a8692a5408", 
                  ask = FALSE, messages = FALSE)[[1]]
  
  # conditions to be satisfied:
  # names should be ok and result is not NULL
  cond <- all(names(res) == c("status", "rank", "division", "scientificname", "commonname", "uid", "genus", "species", "subsp", "modificationdate")) & !is.null(res)
  # if condition is satisfied, print result to file
  if(cond) {
    write.table(x = res, file = "taxid_table.txt", sep = "|", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  } else { print(paste("Check", speclist[i])) }
}

# read-in taxid_table.txt
taxid_df <- read.table(file = "taxid_table.txt", sep = "|", col.names = c("status", "rank", "division", "scientificname", "commonname", "uid", "genus", "species", "subsp", "modificationdate"))
# species with no taxid yet
speclist <- setdiff(unique(df$taxon), unique(taxid_df$scientificname))

# subset df to get representative accession per species that needs taxid
df_sub <- df[df$taxon %in% speclist,]
# remove species that do not have accessions
df_sub <- df_sub[df_sub$genbank_accession != "", ]
# get just one accession per species
df_sub <- df_sub[!duplicated(df_sub$taxon),]
# add empty column for taxids
df_sub$uid <- ""

# use genbanktouid function to get taxids for other species
# 5369.11 eslapsed seconds
for(i in 1:length(df_sub$genbank_accession)) {
  df_sub$uid[i] <- genbank2uid(id = df_sub$genbank_accession[i], 
                               key = "81af93ade140a45a355b621d22a8692a5408")[[1]][1]
}

# check results
df_sub[is.na(df_sub$uid),]
nrow(df_sub[is.na(df_sub$uid),])









##### run codes here:
ptm <- proc.time()
# put function here:
proc.time() - ptm
beepr::beep(5)
