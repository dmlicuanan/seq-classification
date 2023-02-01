# content:
# 1) creates a list of species using a WoRMS dataset
# 2) extracts taxids in MARES database that are absent in WoRMS to create taxa list
# 3) appends taxonomic metadata to WoRMS and MARES taxa lists
# 4) merges WoRMS and MARES metadata to NCBI format

# load packages
library(readr)
library(stringr)
library(data.table)

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

# attach taxonomic metadata to wormstaxlist
dt <- data.table(taxon = wormstaxlist)
# match metadata based on term
dt$taxon_rank <- tolower(d$taxonRank[match(dt$taxon, d$term)])
dt$kingdom <- d$kingdom[match(dt$taxon, d$term)]
dt$phylum <- d$phylum[match(dt$taxon, d$term)]
dt$class <- d$class[match(dt$taxon, d$term)]
dt$order <- d$order[match(dt$taxon, d$term)]
dt$family <- d$family[match(dt$taxon, d$term)]
dt$genus <- d$genus[match(dt$taxon, d$term)]
# add source of metadata
dt$taxon_source <- "worms"
dt$md_source <- "worms"

# keep unique rows only
dt <- unique(dt)
# save RDS to later merge with taxa taken from MARES
# saveRDS(dt, "D:/Documents/NGS/out/Manila_Bay/transients/worms_taxlist_metadata.RDS")





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



# attach taxonomic metadata to marestaxlist
# read-in metadata of wormstaxlist to check columns needed
worms <- readRDS("D:/Documents/NGS/out/Manila_Bay/transients/worms_taxlist_metadata.RDS")
names(worms)
# "taxon"        "taxon_rank"   "kingdom"      "phylum"       "class"        "order"       
# "family"       "genus"        "taxon_source" "md_source"   

# load long taxonomy table (generated from kraken_vis.R) to assign metadata to taxids from mares
long <- readRDS("D:/Documents/NGS/entrez/taxonomy/taxonomy_long.RDS")
# load function for getting metadata (taken from kraken_vis.R)
taxget <- function(ktaxid, rank) {
  linker <- long[parent_taxid == ktaxid, ]$taxid[1]
  taxon <- long[taxid == linker & parent_rank == rank, ]$parent_taxon
  if (length(taxon) == 1) { return(taxon) } else { return(NA) }
}

# create table of taxids from mares not represented in wormstaxlist
mares <- data.table(df[, c("taxid", "taxon")])
rm(df)
# add metadata (all from NCBI)
mares$taxon_rank <- long$parent_rank[match(mares$taxid, long$parent_taxid)]
mares$superkingdom <- sapply(X = mares$taxid, function(x) taxget(x, "superkingdom"))
mares$kingdom <- sapply(X = mares$taxid, function(x) taxget(x, "kingdom"))
mares$phylum <- sapply(X = mares$taxid, function(x) taxget(x, "phylum"))
mares$class <- sapply(X = mares$taxid, function(x) taxget(x, "class"))
mares$order <- sapply(X = mares$taxid, function(x) taxget(x, "order"))
mares$family <- sapply(X = mares$taxid, function(x) taxget(x, "family"))
mares$genus <- sapply(X = mares$taxid, function(x) taxget(x, "genus"))
# add source of metadata
mares$taxon_source <- "mares"
mares$md_source <- "ncbi"
# save RDS to merge with metadata of wormstaxlist
# saveRDS(mares, "D:/Documents/NGS/out/Manila_Bay/transients/mares_taxlist_metadata.RDS")





# prepare worms metadata table (convert to NCBI format)
# load metadata tables
worms <- readRDS("D:/Documents/NGS/out/Manila_Bay/transients/worms_taxlist_metadata.RDS")

# load nodes.RDS and names.RDS to check NCBI format
long <- readRDS("D:/Documents/NGS/entrez/taxonomy/taxonomy_long.RDS")
names <- readRDS("D:/Documents/NGS/entrez/taxonomy/names.RDS")

# taxonomic information from WORMS is inconsistent with NCBI
# for example:
unique(worms$kingdom)
unique(long[parent_rank == "kingdom"]$parent_taxon)
# hence, we use the genus column from WORMS (which was just appended to speciesEpithet) to link each entry in WORMS taxlist to other (higher) taxonomic information

# are all WORMS genera in the long taxonomy table? 15,634 genera are not
gen <- unique(worms$genus)
duds <- gen[!(gen %in% unique(long$parent_taxon))]
# what % species will be "lost" by excluding dud genera?
worms[genus %in% duds] # 47,660 species
(nrow(worms[genus %in% duds]) / nrow(worms))*100 # 21.46189%

# are all genera not in long taxonomy table also absent in names.dmp? 
names[name %in% duds] # 76 genera in names, but not in long
names[(name %in% duds) & (name_class == "synonym")] # 67 genera are synonyms
# create synonym table
st <- names[(name %in% duds) & (name_class == "synonym"), -c(3,4)]
names(st) <- c("taxid", "synonym")
# add column for scientific genera
# subset names to only scientific names
names_sub <- names[name_class == "scientific name"]
st$genus <- names_sub$name[match(st$taxid, names_sub$taxid)]
# what % species will be "lost" if synonym genera are used?
worms[genus %in% setdiff(duds, st$synonym)] # 47521
(nrow(worms[genus %in% setdiff(duds, st$synonym)]) / nrow(worms))*100  # 21.3993%

# add genus_ncbi (as linker to all other taxoxnomic info to worms)
worms$genus_ncbi <- worms$genus
worms[genus %in% st$synonym]$genus_ncbi <- st$genus[match(worms[genus %in% st$synonym]$genus, st$synonym)]
# add taxid of genus_ncbi
worms$genus_taxid <- long$parent_taxid[match(worms$genus_ncbi, long$parent_taxon)]

# subset worms to remove species whose genera are not in long taxonomy table nor have synonyms in names.RDS
worms <- worms[!is.na(genus_taxid)]

# load taxget function
# reduce size of taxonomy table
imp <- unique(long[parent_taxid %in% unique(worms$genus_taxid)]$taxid)
long <- long[taxid %in% imp]

# get unique genus_taxid from worms to reduce running time of taxget
sub <- data.table(genus_taxid = unique(worms$genus_taxid))
# fill-in sub table with taxonomic information
sub$superkingdom2 <- sapply(X = sub$genus_taxid, function(x) taxget(x, "superkingdom"))
sub$kingdom2 <- sapply(X = sub$genus_taxid, function(x) taxget(x, "kingdom"))
sub$phylum2 <- sapply(X = sub$genus_taxid, function(x) taxget(x, "phylum"))
sub$class2 <- sapply(X = sub$genus_taxid, function(x) taxget(x, "class"))
sub$order2 <- sapply(X = sub$genus_taxid, function(x) taxget(x, "order"))
sub$family2 <- sapply(X = sub$genus_taxid, function(x) taxget(x, "family"))

# fill-in worms table with info from sub
worms$superkingdom2 <- sub$superkingdom2[match(worms$genus_taxid, sub$genus_taxid)]
worms$kingdom2 <- sub$kingdom2[match(worms$genus_taxid, sub$genus_taxid)]
worms$phylum2 <- sub$phylum2[match(worms$genus_taxid, sub$genus_taxid)]
worms$class2 <- sub$class2[match(worms$genus_taxid, sub$genus_taxid)]
worms$order2 <- sub$order2[match(worms$genus_taxid, sub$genus_taxid)]
worms$family2 <- sub$family2[match(worms$genus_taxid, sub$genus_taxid)]

# trim worms of columns not needed
worms_ncbi <- worms[, c(12,1,2,13:18,11,9:10)]
# change source of metadata
worms_ncbi$md_source <- "ncbi"

# read-in mares metadata
mares <- readRDS("D:/Documents/NGS/out/Manila_Bay/transients/mares_taxlist_metadata.RDS")

# set preferred variable names
vars <- c("linker_taxid", "taxon", "taxon_rank", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "taxon_source", "md_source")
# check if variable names match
data.frame(vars, names(worms_ncbi), names(mares))
# change variable names of worms and mares
names(worms_ncbi) <- vars
names(mares) <- vars

# bind metadata of wormstaxlist and marestaxlist
l <- list(worms_ncbi, mares)
comp <- rbindlist(l, use.names = TRUE)

# save RDS of compiled taxa from worms and mares, with metadata
saveRDS(comp, "D:/Documents/NGS/out/Manila_Bay/transients/worms_mares_taxlist_metadata.RDS")





# timer
# time process
ptm <- proc.time()
#your function here

proc.time() - ptm