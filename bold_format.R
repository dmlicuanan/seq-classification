# script for 
# 1) formatting BOLD .fastas to kraken format
# 2) deduplicating worms_coi.fasta (which contains both BOLD and NCBI sequences)
# 3) extracting sequences in MARES database not unrepresented in worms_coi.fasta 

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





# formatting BOLD .fastas to kraken format
# set working directory to location of BOLD .fastas
setwd("D:/Documents/NGS/entrez/bold_fastas")
# non-fasta files in folder
list.files()[grep("_bold.fasta", invert = TRUE, list.files())]
# BOLD fastas to process:
files <- list.files(pattern = "_bold.fasta")

# extract ID lines and write to file
# for(i in 1:length(files)) {
#   # read fasta
#   temp <- readFasta(files[i])
#   # write ID to file
#   write(temp$id, file = "bold_id_list", append = TRUE)
# }

# read in file with BOLD ids (453635 rows)
df <- read.table("bold_id_list", sep = "|", fill = TRUE)
# add names
names(df) <- c("processid", "taxon", "markercode", "genbank_accession")





# attempt to attach taxids to BOLD seqs using taxize package

# load required packages
library(taxize)

# species list that need taxids
speclist <- unique(df$taxon)

# 1) generate table for taxids using taxon name; get_uid_ function of taxize
# elspased time: 5.83315 hrs)
# for(i in 1:length(speclist)) {
#   # get taxize::get_uid_ result
#   res <- get_uid_(sci_com = speclist[i], 
#                   key = "81af93ade140a45a355b621d22a8692a5408", 
#                   ask = FALSE, messages = FALSE)[[1]]
#   
#   # conditions to be satisfied:
#   # names should be ok and result is not NULL
#   cond <- all(names(res) == c("status", "rank", "division", "scientificname", "commonname", "uid", "genus", "species", "subsp", "modificationdate")) & !is.null(res)
#   # if condition is satisfied, print result to file
#   if(cond) {
#     write.table(x = res, file = "taxid_table.txt", sep = "|", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
#   } else { print(paste("Check", speclist[i])) }
# }

# read-in taxid_table.txt (21647 rows)
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

# 2) use genbanktouid function to get taxids for other species
# 1.491419 eslapsed hours
for(i in 1:length(df_sub$genbank_accession)) {
  df_sub$uid[i] <- genbank2uid(id = df_sub$genbank_accession[i], 
                               key = "81af93ade140a45a355b621d22a8692a5408")[[1]][1]
}
# write result to table
# write.table(df_sub, file = "taxid_table_genbank2uid.txt", sep = "|", row.names = FALSE, col.names = TRUE, quote = FALSE)

# check results
df_sub <- read.table("taxid_table_genbank2uid.txt", sep = "|", header = TRUE)
# are there blanks in the uid column? no
df_sub[which(df_sub$uid == ""), ]
# are there NAs in the uid column? 39
df_sub[is.na(df_sub$uid),]

# read table of taxids extracted from GenBank files
# codes to generate table in gb_parse.R
# read resulting table
taxid_table_gb <- read.table(file = "D:/Documents/NGS/entrez/wormstaxlist_gb/compiled_gb/taxid_table_genbank.txt", sep = "|", comment.char = "", header = TRUE)

# modify file name to match organism column
taxid_table_gb$file <- gsub("_ncbi.gb", "", taxid_table_gb$file)
taxid_table_gb$file <- gsub("_", " ", taxid_table_gb$file)
# restrict results to organisms not matching with file names
syn <- taxid_table_gb[taxid_table_gb$file != taxid_table_gb$organism, c("organism", "file")]
syn <- syn[!duplicated(syn),]





# checking of taxids from different sources
# clear environment
rm(list = ls())

# re-read relevant tables:
# 1) information of bold taxas
df <- read.table("bold_id_list", sep = "|", fill = TRUE, col.names =  c("processid", "taxon", "markercode", "genbank_accession"))
# taxid tables:
# 2) taxid table from taxize: get_uid_
t1 <- read.table(file = "taxid_table.txt", sep = "|", col.names = c("status", "rank", "division", "scientificname", "commonname", "uid", "genus", "species", "subsp", "modificationdate"))
# 3) taxid table from taxize: genbank2uid
t2 <- read.table("taxid_table_genbank2uid.txt", sep = "|", header = TRUE)
# 4) taxids mined from .gb files
t3 <- read.table(file = "D:/Documents/NGS/entrez/wormstaxlist_gb/compiled_gb/taxid_table_genbank.txt", sep = "|", comment.char = "", header = TRUE)

# fill in bold dataframe with taxids from t1, t2, t3
df$t1 <- t1$uid[match(df$taxon, t1$scientificname)]
df$t2 <- t2$uid[match(df$taxon, t2$taxon)]
df$t3 <- t3$taxid[match(df$taxon, t3$organism)]






# use names.dmp from NCBI taxdump to attach taxids to BOLD seqs
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

# set working directory to location of taxdump
setwd("D:/Documents/NGS/entrez/taxonomy")

# preview files
readLines("readme.txt")
readLines("names.dmp", 10)

# load packages needed
library(data.table)

# read names dump
raw <- fread("names.dmp", sep = "|")
# remove last column
raw <- raw[,-5]

# change variable names
names(raw) <- c("taxid", "name", "name_unique", "name_class")
# variables are:
# tax_id = the id of node associated with this name
# name_txt = name itself
# unique name = the unique variant of this name if name not unique
# name class = (synonym, common name, ...)

# remove tabs from all columns
ndf <- as.data.frame(lapply(raw, trimws))
# remove raw from environment
rm(raw)

# read in bold id list 
bold <- read.table("D:/Documents/NGS/entrez/bold_fastas/bold_id_list", sep = "|", fill = TRUE, col.names =  c("processid", "taxon", "markercode", "genbank_accession"))

# add taxid column
bold$taxid <- ndf$taxid[match(bold$taxon, ndf$name)]

# compare taxids from taxdump with taxids from GenBank files
# read-in table with GenBank taxids
gb <- read.table(file = "D:/Documents/NGS/entrez/wormstaxlist_gb/compiled_gb/taxid_table_genbank.txt", sep = "|", comment.char = "", header = TRUE)
# add column to bold with GenBank taxids
bold$taxid2 <- as.character(gb$taxid[match(bold$taxon, gb$organism)])
# which taxids are not consistent
bold[which(bold$taxid != bold$taxid2),]
# are there taxids in NCBI dump not in GenBank files?
setdiff(unique(bold$taxid2), unique(bold$taxid))

# validate subset to check taxids in https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
set.seed(17)
temp <- unique(bold[,c("taxon", "taxid")])
temp[sample(1:nrow(temp), 20),]

# check number of BOLD sequences with no taxid
notax <- temp[is.na(temp$taxid),]
# number of species with BOLD seq but no taxid = 1075
nrow(notax)
# number of BOLD seqs with no taxid = 4941
nrow(bold[bold$taxon %in% notax$taxon,])
# percentage of BOLD seqs with no taxid = 1.089202%
(nrow(bold[bold$taxon %in% notax$taxon,]) / nrow(bold))*100

# change working directory to where BOLD fastas to be converted are
setwd("D:/Documents/NGS/entrez/bold_fastas")
# set path where converted new fastas (for kraken) will go
path <- "D:/Documents/NGS/entrez/bold_fastas/bold_fastas_kraken/"

# list of bold fastas to be converted to kraken format
files <- list.files(path = "D:/Documents/NGS/entrez/bold_fastas", pattern = "_bold.fasta")

# convert each bold fasta to kraken format (2.481581 hrs elapsed)
for(i in 1:length(files)) {
  # read fasta
  temp <- readFasta(files[i])
  
  # remove gaps from sequences
  temp$seq <- gsub("-", "", temp$seq)
  
  # add species column
  temp$taxon <- sapply(strsplit(temp$id, split = "\\|"), "[[", 2)
  # add taxid
  temp$taxid <- ndf$taxid[match(temp$taxon, ndf$name)]
  
  # remove sequences with no taxid
  temp <- temp[!is.na(temp$taxid),]
  
  # proceed only if temp still has contents
  if(nrow(temp) > 0) {
    # overwrite id column with new id (kraken format)
    temp$id <- paste0("kraken:taxid|", temp$taxid)
    
    # full path of new fasta
    f <- paste0(path, files[i])
    # write fasta
    writeFasta(id = temp$id, seq = temp$seq, file = f)
  }
}

# check that species with no kraken fastas do not have taxids
temp <- setdiff(unique(files), unique(list.files(path = path, pattern = "_bold.fasta")))
# remove file name format
temp <- gsub("_bold.fasta", "", temp)
temp <- sub("_", " ", temp)
# add taxid column
temp <- data.frame(taxon = temp)
temp$taxid <- ndf$taxid[match(temp$taxon, ndf$name)]
# do all species have no associated taxid?
all(is.na(temp$taxid))





# de-duplicate fasta
# set working directory
setwd("D:/Documents/NGS/entrez/kraken_fastas/")

# read in compiled fastas
df <- readFasta("worms_coi.fasta")

# remove duplicates
rdf <- unique(df)
# number of unique species represented = 33156
length(unique(rdf$id))

# write fasta
# writeFasta(id = rdf$id, seq = rdf$seq, file = "worms_coi_deduplicated.fasta")

# compare worms_coi_deduplicated.fasta with MARES pre-compiled database
# read-in pre-compiled database
mares <- readFasta("D:/Documents/NGS/mares/mares_nobar_taxonomy/MARES_NOBAR_BOLD_NCBI_sl_kraken.fasta")

# number of unique species represented in MARES database = 78844
length(unique(mares$id))

# check number of taxids not in worms database = 52621
miss <- gsub("kraken:taxid\\|", "", setdiff(unique(mares$id), unique(rdf$id)))
length(miss)

# check species of taxids in MARES database but not in WORMS database
# load packages needed
library(data.table)

# read names dump
raw <- fread("D:/Documents/NGS/entrez/taxonomy/names.dmp", sep = "|")
# remove last column
raw <- raw[,-5]

# change variable names
names(raw) <- c("taxid", "name", "name_unique", "name_class")
# variables are:
# tax_id = the id of node associated with this name
# name_txt = name itself
# unique name = the unique variant of this name if name not unique
# name class = (synonym, common name, ...)

# remove tabs from all columns
ndf <- as.data.frame(lapply(raw, trimws))
# remove raw from environment
rm(raw)

# dataframe of taxa not in worms_coi database
iden <- ndf[ndf$taxid %in% miss,]
# number of unique taxids in iden = 51631
length(unique(iden$taxid))
# which taxids are listed in the mares database but not in the taxdump
setdiff(miss, unique(iden$taxid))
length(setdiff(miss, unique(iden$taxid)))
# 990 taxids that are in mares but not in the NCBI taxdump nor website






# extracting sequences in MARES database not unrepresented in worms_coi.fasta 
# read in MARES database fasta
mares <- readFasta("D:/Documents/NGS/mares/mares_nobar_taxonomy/MARES_NOBAR_BOLD_NCBI_sl_kraken.fasta")
# add taxid column
mares$taxid <- gsub("kraken:taxid\\|", "", mares$id)

# check mares database
# tabulation of length of id -- all within expected length
table(nchar(unique(mares$id)))
# are there blank sequences?
mares[is.na(mares$seq) | mares$seq == "" | is.na(mares$id),]
# add sequence length column
mares$seqlen <- nchar(mares$seq)
# show sequences with length < 50
mares[mares$seqlen < 50,]

# set directory
setwd("D:/Documents/NGS/entrez/kraken_fastas")

# read in worms fasta database
worms <- readFasta("worms_coi_deduplicated.fasta")
# add taxid column
worms$taxid <- gsub("kraken:taxid\\|", "", worms$id)

# taxids in MARES not in worms fasta
miss <- setdiff(unique(mares$taxid), unique(worms$taxid))

# subset mares so that it contains only taxids not in worms fasta
sub <- mares[mares$taxid %in% miss, ] # 801637
# remove sequences greater than 5000
sub <- sub[sub$seqlen < 2000,]

# write fasta of sequences with taxids in MARES but not in worms fasta
writeFasta(id = sub$id, seq = sub$seq, file = "mares_no_bar_subset.fasta")




##### run codes here:
ptm <- proc.time()
# put function here:

proc.time() - ptm
beepr::beep(5)
