# contents:
# 1. creates (long) taxonomy table from NCBI taxdump so that species taxid can be linked to higher orders of classification
# 2. attempt generating (wide) taxonomy table
# 3. check coverage of taxalists from worms and mares 
# 4. check coverage of COI and 16S databases used for kraken classification
# 5. previews standard output and report files from kraken
# 6. compile results of kraken classficiation per read
# 7. check resolution of kraken classification
# 8. attach other taxonomic information to kraken classifications
# 9. plot NMDS (site ~ taxonomic level)
# 10. check variation between replicates 
# 11. check taxa that are unique, shared, common across replicates
# 12. remove reads present in blanks and extraction negative
# 13. check taxa present in replicates

# load packages
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(vegan)

# load functions
# 1) taxget: add taxa information for kraken classification
# note: taxonomy table long should be loaded
taxget <- function(ktaxid, rank) {
  linker <- long[parent_taxid == ktaxid, ]$taxid[1]
  taxon <- long[taxid == linker & parent_rank == rank, ]$parent_taxon
  if (length(taxon) == 1) { return(taxon) } else { return(NA) }
}

# read-in files for adding other taxonomic information (long taxonomy table)
long <- readRDS("D:/Documents/NGS/entrez/taxonomy/taxonomy_long.RDS")





# taxonomy table (long)
# reformat nodes.dmp to table form (where species/subspecies can be linked to all taxonomic information from rank to superkingdom)
# readme.txt for details:
readLines("D:/Documents/NGS/entrez/taxonomy/readme.txt")

# read nodes.dmp from NCBI taxdump 
nodes <- fread("D:/Documents/NGS/entrez/taxonomy/nodes.dmp", sep = "|", drop = c(7:10,12:14), col.names = c("taxid", "parent", "rank", "emblcode", "divid", "divflag", "hidflag"))
nodes <- as.data.table(lapply(nodes, trimws))
# save nodes RDS
# saveRDS(nodes, "D:/Documents/NGS/entrez/taxonomy/nodes.RDS")

# read names.dmp for assigning taxon names to taxids
names <- fread("D:/Documents/NGS/entrez/taxonomy/names.dmp", sep = "|", drop = 5, col.names = c("taxid", "name", "name_unique", "name_class"))
names <- as.data.table(lapply(names, trimws))
# save RDS
# saveRDS(names, "D:/Documents/NGS/entrez/taxonomy/names.RDS")
# restrict values to scientific names
names <- names[names$name_class == "scientific name", ]

# function that will return parent taxid and parent rank given a taxid 
taxup <- function(taxid) {
  list <- list()
  list$parent_taxid <- nodes$parent[match(taxid, nodes$taxid)]
  list$parent_rank <- nodes$rank[match(list$parent_taxid, nodes$taxid)]
  return(list)
}

# create transition list for reformatting nodes
# elements will be values for each use of taxup
# first element are taxids in nodes.dmp in level of interest:
level <- c("species", "subspecies", "genus", "strain", "no rank")
x <- list(nodes[nodes$rank %in% level,]$taxid) 

# for each taxid, search parent then add to succeeding element in list
# taxup until all taxid values are equal to 1 (corresponds to root)
i <- 1
while( !(all(x[[i]] == "1")) ) {
  # get parent taxid
  x[[i+1]] <- taxup(x[[i]])$parent_taxid
  # get parent rank and add as name
  # names(x[[i+1]]) <- taxup(x[[i]])$parent_rank
  # print progress 
  print(i)
  i <- i + 1
}

# convert list to data.table
dt <- as.data.table(x)
# convert from wide to long format, where 1st column is species/subspecies/etc. taxid
dt <- melt(dt, id.vars = 1, measure.vars = 1:37)
# remove variable column
dt <- dt[ ,-2]
# rename columns
names(dt) <- c("taxid", "parent_taxid")
# remove duplicate columns
dt <- unique(dt)
# add parent rank
dt$parent_rank <- nodes$rank[match(dt$parent_taxid, nodes$taxid)]
# sort data.table by taxid
dt <- dt[order(as.integer(taxid))]

# add taxon names
dt$taxon <- names$name[match(dt$taxid, names$taxid)]
dt$parent_taxon <- names$name[match(dt$parent_taxid, names$taxid)]

# rearrange columns
dt <- dt[ , c(1,4,2,5,3)]
# save as RDS
# saveRDS(dt, "D:/Documents/NGS/entrez/taxonomy/taxonomy_long.RDS")
# save as .txt table
# write.table(dt, file = "D:/Documents/NGS/entrez/taxonomy/taxonomy_long.txt", sep = "\t", row.names = FALSE)





# code below for making taxonomy table wide but:
# not sure if table should still be made wide given that some ranks do not hold unique values (e.g. "clade", "no rank")
# a species can belong to more than one clade
# read RDS
dt <- readRDS("D:/Documents/NGS/entrez/taxonomy/taxonomy_long.RDS")

# convert long to wide: Aggregate function missing, defaulting to 'length'
w <- data.table::dcast(dt, taxid + taxon ~ parent_rank, value.var = "parent_taxon")
table(w$`no rank` > 1)
table(w$clade > 1)
table(w$species > 1)

# restrict table values to those with unique ranks
x <- dt[dt$parent_rank != "no rank" & dt$parent_rank != "clade", ]
# wide format
wide <- data.table::dcast(x, taxid + taxon ~ parent_rank, value.var = c("parent_taxon", "parent_taxid"))
# arrange by taxid
wide <- wide[order(as.integer(taxid))]
# save .txt table
# write.table(wide, file = "D:/Documents/NGS/entrez/taxonomy/taxonomy_wide.txt", sep = "\t", row.names = FALSE, quote = TRUE)
# save RDS
# saveRDS(wide, "D:/Documents/NGS/entrez/taxonomy/taxonomy_wide.RDS")





# check coverage of taxalists from worms and mares
# load taxalists with metadata
raw <- readRDS("D:/Documents/NGS/out/Manila_Bay/transients/worms_mares_taxlist_metadata.RDS")
# remove entries with no associated taxon (some taxids from mares not in taxdumps)
md <- raw[!is.na(taxon)]

# make taxon_source a factor
md$taxon_source <- factor(md$taxon_source, levels = c("worms", "mares"))
# remove linker_taxid and md_source columns and remove duplicate rows
md <- md[ , -c("linker_taxid", "md_source")]
md <- md[!duplicated(md[,-c("taxon_source")])]
# check which are duplicates: 
# all are okay except taxa from worms whose metadata are derived from genus
# Gaimardia, Kali and Squilla are plant AND animal genera, so we remove plant entries
md[taxon %in% md[duplicated(taxon)]$taxon][order(taxon)]
rem <- md[taxon %in% md[duplicated(taxon)]$taxon & taxon_source == "worms"]
# remove entries in rem from md
md <- md[!(taxon %in% rem$taxon & taxon_source == "worms")]

# number of taxa in worms and mares taxlists combined
length(unique(md$taxon))
# number of echinoderm species represented
ech <- md[phylum == "Echinodermata" & taxon_rank == "species"]
length(unique(ech$taxon))
table(ech$class)

# restrict to only Eukaryotic phyla
sub <- md[superkingdom == "Eukaryota" & !is.na(phylum) & taxon_rank == "species"]
# determine order by which bars will be arranged
temp <- sub[, .N, by = .(phylum)][order(N)]
# distribution of taxa in taxlists across phyla
data <- sub[, .N, by = .(phylum, kingdom, taxon_source)][order(N)]
# make kingdom and phylum factors
data$phylum <- factor(data$phylum, levels = temp$phylum)
data$kingdom <- gsub("Viridiplantae", "Plants", data$kingdom)
data$kingdom <- factor(data$kingdom, levels = c("Metazoa", "Plants", "Fungi"))
# grouped, stacked barplot
ggplot(data, aes(fill = taxon_source, y = N, x = phylum)) + 
  geom_bar(position= "stack", stat = "identity", width = 0.5) + 
  scale_y_continuous(breaks = seq(0,67000,6700)) +
  ylab("Number species taxids") +
  xlab("Phylum") +
  labs(fill = "Taxon source") +
  theme_minimal() + coord_flip() +
  facet_grid(kingdom ~ ., scales = "free_y", space = "free_y", switch = NULL)

# restrict to non-Eukaryotic taxa
sub <- md[superkingdom %in% c("Archaea",  "Bacteria", "Viruses") & !is.na(phylum)]
# determine order by which bars will be arranged
temp <- sub[, .N, by = .(phylum)][order(N)]
# distribution of taxa in taxlists across phyla
data <- sub[, .N, by = .(phylum, superkingdom, taxon_source)][order(N)]
# make phylum factor
data$phylum <- factor(data$phylum, levels = temp$phylum)
# grouped, stacked barplot
ggplot(data, aes(fill = taxon_source, y = N, x = phylum)) + 
  geom_bar(position= "stack", stat = "identity", width = 0.5) + 
  scale_y_continuous(breaks = seq(0,350,35)) +
  ylab("Number of taxids") +
  xlab("Phylum") +
  labs(fill = "Taxon source") +
  theme_minimal() + coord_flip() +
  facet_grid(superkingdom ~ ., scales = "free_y", space = "free_y", switch = NULL)





# check coverage of COI and 16S databases used for kraken classification
# load databases needed 
db1 <- readLines("D:/Documents/NGS/entrez/kraken_fastas/worms_mares_coi.fasta")
db2 <- readLines("D:/Documents/NGS/entrez/kraken_fastas/worms_mares_16S_deduplicated.fasta")
# 2 lines per fasta, so number of fastas per database is:
sapply(X = list(db1, db2), function(x) length(x)/2)

# restrict to lines with "taxid" string and remove ">kraken:taxid|" string
db1 <- gsub("^>kraken:taxid\\|", "", db1[grepl(":taxid", db1)])
db2 <- gsub("^>kraken:taxid\\|", "", db2[grepl(":taxid", db2)])

# tabulate number of fastas per taxid
db1 <- as.data.table(table(db1))
db2 <- as.data.table(table(db2))
# number of taxids per database
sapply(X = list(db1, db2), nrow)

# read-in files for adding other taxonomic information (long taxonomy table)
long <- readRDS("D:/Documents/NGS/entrez/taxonomy/taxonomy_long.RDS")
# read-in original nodes and names
# nodes <- readRDS("D:/Documents/NGS/entrez/taxonomy/nodes.RDS")
# names <- readRDS("D:/Documents/NGS/entrez/taxonomy/names.RDS")
# add rank of taxon 
long$taxon_rank <- long$parent_rank[match(long$taxid, long$parent_taxid)]
# number of echinoderm species in listed in NCBI: 5281
long[parent_taxon == "Echinodermata" & taxon_rank == "species"]

# merge fasta counts in one data.table
# create taxid column based on taxids of coi and 16S fastas
dt <- data.table(taxid = unique(c(db1$db1, db2$db2)))
# fill in coi and sxs columns with fasta counts
dt$coi_N <- db1$N[match(dt$taxid, db1$db1)]
dt$sxs_N <- db2$N[match(dt$taxid, db2$db2)]
# add rank information
dt$rank <- long$parent_rank[match(dt$taxid, long$parent_taxid)]
# remove objects not needed
rm(db1, db2)
# reduce size of taxonomy table since it takes up too much space
long <- long[(taxid %in% dt$taxid) | (parent_taxid %in% dt$taxid)]

# number of taxids unique to each database
nrow(dt[coi_N > 0 & is.na(sxs_N)]) # 82012
nrow(dt[sxs_N > 0 & is.na(coi_N)]) # 16519
# subset by database
temp <- dt[coi_N > 0]
temp <- dt[sxs_N > 0]
# taxonomic levels of taxids
table(temp$rank)
pie(table(temp$rank))
pie(table(dt$rank))

# check taxa in databases with no ranks: 1027 taxids
temp <- dt[is.na(rank),]
# taxids not in long taxonomy table but present in nodes.dmp:
# 35 are absent because rank is more specific than species (e.g. varietas, forma, isolate, biotype)
nodes[nodes$taxid %in% temp$taxid, ]
# other 992 remaining taxa with no rank: just absent in taxdump files
dt[is.na(rank) & !(taxid %in% nodes$taxid),]
dt[is.na(rank) & !(taxid %in% names$taxid),]

# add superkingdom, kingdom, phylum, and class to dt
dt$superkingdom <- sapply(X = dt$taxid, function(x) taxget(x, "superkingdom"))
dt$kingdom <- sapply(X = dt$taxid, function(x) taxget(x, "kingdom"))
dt$phylum <- sapply(X = dt$taxid, function(x) taxget(x, "phylum"))
dt$class <- sapply(X = dt$taxid, function(x) taxget(x, "class"))
# save RDS
# saveRDS(dt, "D:/Documents/NGS/out/Manila_Bay/transients/kraken_dbs_info.RDS")

# restrict to Eukaryote species in COI database and remove NA phyla
sub <- dt[(superkingdom == "Eukaryota") & !(is.na(phylum)) & (rank == "species") & (coi_N > 0)]
# count number of species taxids per phylum
sub <- sub[, .N, by = .(phylum, kingdom)][order(N)]
# make kingdom and phylum factors
sub$phylum <- factor(sub$phylum, levels = sub$phylum)
sub$kingdom <- factor(sub$kingdom, levels = c("Metazoa", "Viridiplantae", "Fungi"))
# distribution of fastas across phyla
ggplot(sub, aes(x = phylum, y = N, color = kingdom)) +
  geom_segment(aes(x = phylum, xend = phylum, y = 0, yend = N), color = "black") +
  geom_point(size = 2.5) +
  scale_y_continuous(breaks = seq(0,35000,3500)) +
  xlab("Phylum") + ylab("No. of species represented in COI database") +
  labs(color = "Kingdom") +
  facet_grid(kingdom ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_minimal() +
  theme(strip.text = element_blank()) +
  coord_flip() 

# restrict to Eukaryote species in 16S database and remove NA phyla
sub <- dt[(superkingdom == "Eukaryota") & !(is.na(phylum)) & (rank == "species") & (sxs_N > 0)]
# count number of species taxids per phylum
sub <- sub[, .N, by = .(phylum, kingdom)][order(N)]
# make kingdom and phylum factors
sub$phylum <- factor(sub$phylum, levels = sub$phylum)
sub$kingdom <- factor(sub$kingdom, levels = c("Metazoa", "Viridiplantae", "Fungi"))
# distribution of fastas across phyla
ggplot(sub, aes(x = phylum, y = N, color = kingdom)) +
  geom_segment(aes(x = phylum, xend = phylum, y = 0, yend = N), color = "black") +
  geom_point(size = 2.5) +
  scale_y_continuous(breaks = seq(0,3500,350)) +
  xlab("Phylum") + ylab("No. of species represented in 16S database") +
  labs(color = "Kingdom") +
  facet_grid(kingdom ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_minimal() +
  theme(strip.text = element_blank()) +
  coord_flip() 
# references:
# https://www.data-to-viz.com/caveat/pie.html
# https://stackoverflow.com/questions/71263664/plot-values-with-a-certain-order-based-on-another-column-in-ggplot

# echinoderm counts
# number of taxids (species, genus, family, subspecies) under Echinodermata: 1681
# number of species taxids
dt[phylum == "Echinodermata" & rank == "species"] # 1457
# subset based on database
temp <- dt[phylum == "Echinodermata" & rank == "species" & coi_N > 0] # 1456
temp <- dt[phylum == "Echinodermata" & rank == "species" & sxs_N > 0] # 42
# number of species fastas
sum(temp$coi_N)
sum(temp$sxs_N)
# distribution of classes
table(temp$class)
# without restriction to rank
temp <- dt[phylum == "Echinodermata" & coi_N > 0]
sum(temp$coi_N)
temp <- dt[phylum == "Echinodermata" & sxs_N > 0]
sum(temp$sxs_N)





# preview kraken standard output and reports
# location of kraken files
setwd("D:/Documents/NGS/out/Manila_Bay/kraken/")

# files that underwent merging then trimming
list.files(pattern = "merged_trimmed") # 136
list.files(pattern = "_merged_trimmed.kraken$") # 34
list.files(pattern = "_merged_trimmed.krakenReport") # 34
list.files(pattern = "unmerged_trimmed.kraken$") # 34
list.files(pattern = "unmerged_trimmed.krakenReport") # 34

# preview kraken standard output and report files
# kraken standard output and report files
f1 <- list.files(pattern = ".kraken$")
f2 <- list.files(pattern = ".krakenReport$")
# index of file to be read
i <- 1

# preview file
readLines(f1[i], 10)
# .kraken files have as many lines are there are reads
# classified reads start with C
# unclassified reads start with U
readLines(f2[i], 10)

# load file to see format
d1 <- fread(f1[i], sep = "\t", col.names = c("code", "seqid", "taxon", "rdlen", "lcamap"))
d2 <- fread(f2[i], sep = "\t", col.names = c("percroot", "ctroot", "ctdirect", "rank", "taxid", "sciname"))





# compile results per read
# location of kraken files
setwd("D:/Documents/NGS/out/Manila_Bay/kraken/")

# list of kraken standard output to be compiled
# includes merged and unmerged
files <- list.files(pattern = "merged_trimmed.kraken$")

# read-in files for adding other taxonomic information (long taxonomy table)
long <- readRDS("D:/Documents/NGS/entrez/taxonomy/taxonomy_long.RDS")

# create empty list
df <- vector(mode = "list", length = length(files))
# read each standard output file to fill-in list
for(i in 1:length(files)) {
  # basic output
  temp <- fread(files[i], sep = "\t", col.names = c("code", "seqid", "taxon", "rdlen", "lcamap"))
  
  # include source file
  temp$file <- files[i]
  
  # tag taxon with site name
  temp$site <- sub("-.*", "", temp$file)
  # tag with primer used
  temp$primer <- ifelse(grepl("-UM_", temp$file), "UM",
                       ifelse(grepl("-16S-", temp$file), "16S", 
                              ifelse(grepl("-OPH_", temp$file), "OPH", NA)))
  # merged or unmerged?
  temp$merged <- ifelse(grepl("_merged", temp$file), "merged", 
                        ifelse(grepl("_unmerged", temp$file), "unmerged", NA))
  # sample, blank, or extraction negative
  temp$type <- ifelse(grepl("^b", temp$file), "blank",
                      ifelse(grepl("^EN", temp$file), "negative", "sample"))
  
  # extract taxid
  temp$taxid <- gsub(".*\\(taxid (.+)\\)", "\\1", temp$taxon)  
  # add temp to list
  df[[i]] <- temp
}

# collapse all kraken stdout in one list
df <- rbindlist(df)

# taxids that are not in taxonomy table: unclassified (taxid 0) 
unique(df$taxid)[!(unique(df$taxid) %in% unique(long$parent_taxid))]

# add rank from taxonomy table
df$taxon_rank <- long$parent_rank[match(df$taxid, long$parent_taxid)]
# which taxa have no rank from taxonomy table
df[is.na(df$taxon_rank), ]
# mark unclassifieds
df[df$taxid == "0", ]$taxon_rank <- "unclassified"

# set order of ranks
ranks <- c("superkingdom", "kingdom", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "cohort", "superorder", "order", "parvorder", "suborder", "infraorder", "superfamily", "family", "subfamily", "tribe", "genus", "subgenus", "species", "subspecies", "strain", "clade", "no rank", "unclassified")
# convert taxon_rank as factor
df$taxon_rank <- factor(df$taxon_rank, levels = ranks)

# add name (from NCBI) 
df$ncbi_taxon <- long$parent_taxon[match(df$taxid, long$parent_taxid)]
# compare with taxon name from kraken
temp <- df[,c("taxon", "ncbi_taxon")]
# remove text between parenthesis
temp$taxon <- gsub("\\s\\(taxid.*\\)", "", temp$taxon)
# non-matching names: all are synonyms
unique(temp[(temp$taxon != temp$ncbi_taxon) | is.na(temp$ncbi_taxon),])

# save as RDS (read per row)
# saveRDS(df, "D:/Documents/NGS/out/Manila_Bay/transients/compiled_kraken_stdout_per_read.RDS")

# read RDS
df <- readRDS("D:/Documents/NGS/out/Manila_Bay/transients/compiled_kraken_stdout_per_read.RDS")

###### convert from per read view to per taxon (per site) view
# df <- df[, .N, by = .(taxon, code, file, site, primer, merged, type, taxid)]

# simplify rank
# set order of ranks
ranks <- c("superkingdom", "kingdom", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "cohort", "superorder", "order", "parvorder", "suborder", "infraorder", "superfamily", "family", "subfamily", "tribe", "genus", "subgenus", "species", "subspecies", "strain", "clade", "no rank", "unclassified")
# simplified ranks
simp <- c("superkingdom", "kingdom", "phylum", "phylum", "class", "class", "class", "class", "class", "order", "order", "order", "order", "order", "family", "family", "family", "family", "genus", "genus", "species", "species", "species", "clade", "no rank", "unclassified")
# create data frame for matching original ranks with simplified ranks
r <- data.frame(ranks, simp)

# add simplified rank column
df$taxon_rank_simp <- r$simp[match(df$taxon_rank, r$ranks)]
# convert as factor
df$taxon_rank_simp <- factor(df$taxon_rank_simp, levels = unique(simp))

# convert type as factor
df$type <- factor(df$type, levels = c("sample", "blank", "negative"))





# resolution of kraken classification
# get number of reads per rank
data <- df[, .N, by = .(taxon_rank_simp, type)][order(-taxon_rank_simp)]
# plot colors
ncolor <- length(unique(data$taxon_rank_simp))
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
# percent stacked barplot
p <- ggplot(data, aes(fill = taxon_rank_simp, y = N, x = type)) + 
  geom_bar(position= "fill", stat = "identity", width = 0.5) + 
  ylab("Proportion of reads classified to rank") +
  labs(fill = "rank") +
  scale_fill_manual(values = getPalette(ncolor)) +
  theme_minimal()
# stacked barpchart
p <- ggplot(data, aes(fill = taxon_rank_simp, y = N, x = type)) + 
  geom_bar(position= "stack", stat = "identity", width = 0.5) + 
  ylab("Number of reads classified to rank") +
  labs(fill = "rank") + 
  scale_fill_manual(values = getPalette(ncolor)) +
  theme_minimal()
# plot
ggplotly(p)

# data summary: check rank composition of reads per primer
temp <- split(data, data$primer)
# 16S and OPH
temp[[1]]$perc <- (temp[[1]]$N/sum(temp[[1]]$N))*100
temp[[2]]$perc <- (temp[[2]]$N/sum(temp[[2]]$N))*100
# UM
um <- split(temp[[3]], temp[[3]]$type)
for(i in 1:3) {
  um[[i]]$perc <- (um[[i]]$N/sum(um[[i]]$N))*100
}

# overall number of taxa per rank
# number of sites per primer
temp <- unique(df[,.(site, primer, type)])
temp[, .(sites = .N), by = .(primer, type)]

# list of taxa per per primer and type
data <- unique(df[,.(taxon, primer, taxon_rank_simp, type)])
# number of unique taxa per rank, primer, and type
data <- data[, .N, by = .(taxon_rank_simp, type, primer)][order(-taxon_rank_simp)]

# selected color palette
pal <- "Set1"

# subset data (samples, excluding blanks and negatives)
temp <- data[type == "sample" & taxon_rank_simp != "unclassified",]
# number of taxa vs. taxonomic rank
ggplot(temp, aes(x = taxon_rank_simp, y = N, group = primer, color = primer)) + 
  geom_line() + geom_point() +
  ylab("Number of taxa") + xlab("Taxonomic rank") +
  labs(color = "Primer") +
  scale_colour_brewer(palette = pal) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  ggtitle("Samples only")

# subset data (blanks and negatives)
temp <- data[(type == "blank" | type == "negative") & taxon_rank_simp != "unclassified",]
# get color palette of primers
display.brewer.pal(3, pal)
col <- brewer.pal(n = 3, pal)[3]
# plot
ggplot(temp, aes(x = taxon_rank_simp, y = N, group = type))  + 
  geom_line(aes(linetype = type), color = col) + 
  geom_point(color = col) +
  ylab("Number of taxa") + xlab("Taxonomic rank") +
  labs(color = "Primer") +
  scale_colour_manual(col) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  ggtitle("Blanks and negative")





# add other taxonomic metadata for kraken classifications
# get list of unique taxids that taxonomic info are needed for
list <- data.table(taxid = unique(df$taxid))
# retrieve taxonomic information of each taxid
list$superkingdom <- sapply(X = list$taxid, function(x) taxget(x, "superkingdom"))
list$kingdom <- sapply(X = list$taxid, function(x) taxget(x, "kingdom"))
list$phylum <- sapply(X = list$taxid, function(x) taxget(x, "phylum"))
list$class <- sapply(X = list$taxid, function(x) taxget(x, "class"))
list$order <- sapply(X = list$taxid, function(x) taxget(x, "order"))
list$family <- sapply(X = list$taxid, function(x) taxget(x, "family"))
list$genus <- sapply(X = list$taxid, function(x) taxget(x, "genus"))

# append taxonomic information to df
df$superkingdom <- list$superkingdom[match(df$taxid, list$taxid)]
df$kingdom <- list$kingdom[match(df$taxid, list$taxid)]
df$phylum <- list$phylum[match(df$taxid, list$taxid)]
df$class <- list$class[match(df$taxid, list$taxid)]
df$order <- list$order[match(df$taxid, list$taxid)]
df$family <- list$family[match(df$taxid, list$taxid)]
df$genus <- list$genus[match(df$taxid, list$taxid)]

# save as RDS (read per row)
# saveRDS(df, "D:/Documents/NGS/out/Manila_Bay/transients/compiled_kraken_stdout_per_read.RDS")





# ordination 
# read-in data
df <- readRDS("D:/Documents/NGS/out/Manila_Bay/transients/compiled_kraken_stdout_per_read.RDS")

# subset data
sub <- df[primer == "UM"]
sub <- df

# summarize number of reads by site
temp <- sub[, .(reads = .N), by = .(site, taxid)]
temp <- sub[, .(reads = .N), by = .(site, phylum)]
temp <- sub[, .(reads = .N), by = .(site, class)]
temp <- sub[, .(reads = .N), by = .(site, family)]
temp <- sub[, .(reads = .N), by = .(site, genus)]

# make columns taxa, and rows the sites
dt <- dcast(temp, site ~ taxid, value.var = "reads")
dt <- dcast(temp, site ~ phylum, value.var = "reads")
dt <- dcast(temp, site ~ class, value.var = "reads")
dt <- dcast(temp, site ~ family, value.var = "reads")
dt <- dcast(temp, site ~ genus, value.var = "reads")
# replace all NAs with 0
dt[is.na(dt)] <- 0

# create presence-absense matrix  
mat <- as.matrix(dt[,-1], rownames = dt$site)
# replace all non-zeros with 1
mat[mat > 0] <- 1

# compute jaccard dissimilarity indices
distmat <- vegdist(mat, method = "jaccard")

# hierarchical clustering based on jaccard
plot(hclust(distmat))

# run NMDS
m <- metaMDS(distmat, display = "sites")
# base plots
plot(m, type = "t")
ordihull(m, rownames(m$points))

# ggplots
# get nmds scores
c <- as.data.frame(scores(m)) 
# add site names
c$site <- rownames(c) 
# add group
c$group <- "negative"
c[grepl("^1|^2|^8|^9", c$site),]$group <- "south"
c[grepl("^3|^4|^5|^6|^7", c$site),]$group <- "north"
# extract site number
c$siteno <- sub("(^[0-9]+).*$", "\\1", c$site)
c$siteno[grepl("^[[:alpha:]]", c$siteno)] <- "negative"

# NMDS based on presence-absence of taxids across sites, regardless of primer used
ggplot(c,aes(x = NMDS1, y = NMDS2, colour = siteno, shape = group)) + 
  geom_point(size=3) +
  geom_text(aes(label = site), colour = "black", size = 3, vjust = -1, hjust = 1) +
  labs(color = "Site", shape = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  coord_equal() +
  scale_color_brewer(palette="Paired")





# group by primer
# summarize number of reads by site and primer
temp <- df[, .(reads = .N), by = .(site, primer, taxid)]
temp <- df[, .(reads = .N), by = .(site, primer, phylum)]
temp <- df[, .(reads = .N), by = .(site, primer, class)]
temp <- df[, .(reads = .N), by = .(site, primer, family)]
temp <- df[, .(reads = .N), by = .(site, primer, genus)]

# paste site and primer
temp$id <- paste0(temp$site, "_", temp$primer)

# make columns taxa, and rows the sites
dt <- dcast(temp, id ~ taxid, value.var = "reads")
dt <- dcast(temp, id ~ phylum, value.var = "reads")
dt <- dcast(temp, id ~ class, value.var = "reads")
dt <- dcast(temp, id ~ family, value.var = "reads")
dt <- dcast(temp, id ~ genus, value.var = "reads")

# replace all NAs with 0
dt[is.na(dt)] <- 0

# create presence-absense matrix  
mat <- as.matrix(dt[,-1], rownames = dt$id)
# replace all non-zeros with 1
mat[mat > 0] <- 1

# compute jaccard dissimilarity indices
distmat <- vegdist(mat, method = "jaccard")

# run NMDS
m <- metaMDS(distmat, display = "sites")
# base plots
plot(m, type = "t")

# ggplots
# get nmds scores
c <- as.data.frame(scores(m)) 
# add id names
c$id <- rownames(c) 
# extract primer and site
c$primer <- gsub(".*_", "", c$id)
c$site <- gsub("_.*", "", c$id)
# extract site number
c$siteno <- sub("(^[0-9]+).*$", "\\1", c$site)
c$siteno[grepl("^[[:alpha:]]", c$siteno)] <- "negative"

# NMDS based on presence-absence of taxids across sites, grouped by primer
ggplot(c,aes(x = NMDS1, y = NMDS2, colour = siteno, shape = primer)) + 
  geom_point(size=3) +
  geom_text(aes(label = site), colour = "black", size = 2.5, vjust = -1, hjust = 1) +
  labs(shape = "Primer", colour = "Site") +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  coord_equal() + 
  scale_color_brewer(palette="Paired")





# check variation between replicates 
# read-in data
df <- readRDS("D:/Documents/NGS/out/Manila_Bay/transients/compiled_kraken_stdout_per_read.RDS")
# restrict entries to those with replicates (samples sequenced using uniminibarcode primer)
# remove site 4 since it has no replicate
# remove unclassified reads
sub <- df[primer == "UM" & type == "sample" & site != "4A" & code == "C"]
# summarize number of reads by site, etc (so that there is one row for each taxon in a replicate)
sub <- sub[, .(reads = .N), by = .(site, taxid, taxon, taxon_rank_simp, superkingdom, kingdom, phylum, class, order, family, genus)]
# add site number
sub$siteno <- gsub("\\D+", "", sub$site)
unique(sub[,c("site", "siteno")])

# create empty list for counts
names <- c("site", "total", "mean", "sem", "common", "shared", "unique")
l <- sapply(names, function(x) NULL)
# create empty data.frame where common, shared, unique taxa will be saved
tl <- data.frame()

# add names of sites with replicates
l$site <- unique(sub$siteno)
# fill-in counts table
for(i in 1:length(l$site)) {
  # create subset for site
  temp <- sub[siteno == l$site[i]]
  # set taxonomic rank to be used for making counts 
  temp$x <- temp$phylum
  
  # remove NAs 
  temp <- temp[!is.na(x)]
  # get unique values of taxon per replicate
  temp <- unique(temp[, c("siteno", "site", "x")])
   
  # subset based on replicate
  reps <- split(temp, temp$site)
  
  # get total richness in site
  l$total[i] <- length(unique(temp$x))
  # get mean richness across replicates and standard error of the mean
  l$mean[i] <- round(mean(temp[, .N, by = .(site)]$N), 1)
  l$sem[i] <- round(sd(temp[, .N, by = .(site)]$N)/sqrt(length(temp[, .N, by = .(site)]$N)), 1)
  
  # there are two replicates
  if(length(reps) == 3) {
    a <- reps[[1]]$x
    b <- reps[[2]]$x
    c <- reps[[3]]$x
    
    # taxa common across all replicates
    common <- Reduce(intersect, list(a, b, c))
    l$common[i] <- length(common)
    
    # taxa shared by at most two replicates
    atleast <- Reduce(union, list(intersect(a, b), intersect(a, c), intersect(b,c)))
    shared <- setdiff(atleast, common)
    l$shared[i] <- length(shared)
    
    # unique 
    unique <- Reduce(union, list(setdiff(a, union(b,c)), setdiff(b, union(a,c)), setdiff(c, union(a,b))))
    l$unique[i] <- length(unique)
    
    # compile taxa in list
    t1 <- data.frame(taxon = common, category = "common")
    t2 <- data.frame(taxon = shared, category = "shared")
    t3 <- data.frame(taxon = unique, category = "unique")
    t <- rbind(t1, t2, t3)
  }
  
  if(length(reps) == 2) {
    a <- reps[[1]]$x
    b <- reps[[2]]$x
    
    # taxa common across all replicates
    common <- intersect(a,b)
    l$common[i] <- length(common)
    
    # NA since there are only two sets
    l$shared[i] <- NA
    
    # unique 
    unique <- union(setdiff(a, b), setdiff(b, a))
    l$unique[i] <- length(unique)
    
    # compile taxa in list
    t1 <- data.frame(taxon = common, category = "common")
    t3 <- data.frame(taxon = unique, category = "unique")
    t <- rbind(t1, t3)
  }
  
  # add site number to taxa list
  t$site <- l$site[i]
  # compile list of taxa and if they are common, shared or unique
  tl <- rbind(tl, t)
}

# common, shared, and unique taxa as data frame
d <- data.frame(l)

# melt data frame for plotting
dt <- melt(data.table(d), id.vars = "site", measure.vars = c("common", "shared", "unique"))
ggplot(dt, aes(fill = variable, y = value, x = site)) + 
  geom_bar(position= "stack", stat = "identity", width = 0.5) + 
  labs(fill = "", x = "Site", y = "Number of taxa") + 
  scale_x_discrete(limits = rev) + 
  scale_fill_brewer(palette="Set3", labels = c("common across all replicates", "shared by at most 2 replicates", "unique to replicate")) +
  geom_text(aes(label = value), colour = "black", position = position_stack(vjust = 0.5), size = 3) +
  theme_minimal() + coord_flip() 

# check taxa that are common, shared, unique
tl <- data.table(tl)
# get unique values per category
unique <- unique(tl[category == "unique"]$taxon)
common <- unique(tl[category == "common"]$taxon)
shared <- unique(tl[category == "shared"]$taxon)
# get unique phyla
all <- unique(c(unique, shared, common))
# convert phyla to factor
all <- factor(all, levels = all)

# get presence-absence of phyla across categories
pa <- data.table(phylum = all)
pa$unique <- pa$phylum %in% unique
pa$shared <- pa$phylum %in% shared
pa$common <- pa$phylum %in% common
# melt data table
pa <- melt(pa, id.vars = 1, measure.vars = 2:4)
names(pa) <- c("phylum", "category", "presence")
# convert presence as factor
pa$presence <- factor(pa$presence, levels = c(TRUE, FALSE))

# create geom_tile (presence-absence matrix)
ggplot(pa, aes(x = category, y = phylum, fill = presence)) +
  geom_tile(colour = "black") + 
  scale_fill_manual(values = c("olivedrab1", "snow2"), labels = c("present", "absent")) +
  scale_y_discrete(limits = rev) + 
  scale_x_discrete(position = "top", labels = c("Unique", "Shared", "Common")) +
  labs(x = "", y = "Phylum", fill = "") + 
  coord_fixed(ratio = 1/5) +
  theme_minimal() + 
  theme(axis.ticks = element_blank())

# edit tl table so that "common" and "shared" are lumped together
# shared will now mean taxa that are shared by at least two replicates
tl$category <- gsub("common", "shared", tl$category)

# add "absent" category
# list of all phyla
phy <- unique(tl$taxon)
# create df of taxa that are absent in every site
emp <- data.frame()
for(i in 1:length(unique(tl$site))) {
  temp <- tl[site == unique(tl$site)[i]]
  emp <- rbind(emp, data.frame(taxon = setdiff(phy, temp$taxon),
                               category = "absent",
                               site = unique(tl$site)[i]))
}
# complete tl with "absent" category
pd <- rbind(tl, emp)
# convert category as factor
pd$category <- factor(pd$category, levels = c("shared", "unique", "absent"))

# original idea for plot:
ggplot(tl, aes(y = taxon, x = site, fill = category)) +
  geom_tile() +
  facet_grid(~category) 
# combined fill
ggplot(pd, aes(y = taxon, x = site, fill = category)) +
  geom_tile(colour = "black") +
  scale_fill_manual(values = c("darkolivegreen2", "firebrick1", "grey90"),
                    labels = c("Shared by at least two replicates",
                               "Unique to one replicate in site",
                               "Undetected from site")) +
  scale_y_discrete(limits = rev) +
  labs(x = "Site", y = "Phylum", fill = "") +
  coord_fixed(ratio = 1/2) +
  theme_minimal() +
  theme(axis.ticks = element_blank())





# remove reads present in blanks and extraction negative
# read-in data
df <- readRDS("D:/Documents/NGS/out/Manila_Bay/transients/compiled_kraken_stdout_per_read.RDS")

# subset df (remove unclassified reads) and restrict to COI samples since only those have blanks and negatives
sub <- df[taxid != "0" & primer == "UM"]

# get number of reads per taxon per replicate (samples)
samp <- sub[type == "sample", .(reads = .N), by = .(site, primer, taxid, taxon, taxon_rank_simp, superkingdom, kingdom, phylum, class, order, family, genus)]
# add corresponding blank per sample
samp$blank <- ifelse(grepl("^1|^2", samp$site), "b1B", 
                     ifelse(grepl("^8|^9", samp$site), "b8A", "b5A"))
# check assigned blanks
unique(samp[,c("site", "blank")])

# get number of reads per taxon per replicate (extraction negative)
neg <- sub[type == "negative", .(reads = .N), by = .(site, taxid)]
# get number of reads per taxon per replicate (blanks)
blanks <- sub[type == "blank", .(reads = .N), by = .(site, taxid)]

# get number of extraction negative reads that will be subtracted
samp$reads_neg <- neg$reads[match(samp$taxid, neg$taxid)]
# which taxa will be removed (those with reads less than or equal to zero)
samp[reads - reads_neg <= 0]

# list of replicates
reps <- unique(samp$site)
# add empty column
samp$reads_blank <- 0
# fill-in number of reads (from blanks) that will be subtracted
for(i in 1:length(reps)) {
  # subset per replicate
  s <- samp[site == reps[i]]
  # subset blanks to corresponding blank of replicate
  b <- blanks[site == unique(s$blank)]
  
  # match number of blank reads
  s$reads_blank <- b$reads[match(s$taxid, b$taxid)]
  samp[samp$site == reps[i]]$reads_blank <- s$reads_blank
}

# replace all NAs with 0
samp$reads_neg[is.na(samp$reads_neg)] <- 0
samp$reads_blank[is.na(samp$reads_blank)] <- 0
# get final reads excluding blanks and negative
samp$reads_fin <- samp$reads - samp$reads_neg - samp$reads_blank

# check number of taxids that will be removed
a <- samp[reads_blank > 0 | reads_neg > 0, c("site", "taxid", "reads", "reads_neg", "reads_blank", "reads_fin")]
a 
table(a$site) # taxa per replicate
table(a[reads_fin > 0]$site) # taxa retained per replicate
table(a[reads_fin <= 0]$site) # taxa to be removed per replicate

# format reads to be retained for merging with reads from other primers
temp <- samp[reads_fin > 0, -(13:16)]

# get reads from other primers
sub <- df[taxid != "0" & primer != "UM"]
# get number of reads per taxon per replicate (samples)
oth <- sub[type == "sample", .(reads_fin = .N), by = .(site, primer, taxid, taxon, taxon_rank_simp, superkingdom, kingdom, phylum, class, order, family, genus)]

# bind remaining reads from UM with untouched reads from 16S and OPH
data <- rbind(oth, temp)
# consider only taxa identified until the class level
# data <- data[taxon_rank_simp == "class"]

# get number of taxa per phylum
pd <- data[, .N, by = .(site, phylum, primer)]
# get the top most common phyla per primer
a <- pd[primer == "16S" & !is.na(phylum), .N, by = .(phylum)][order(-N)]$phylum[1:10]
b <- pd[primer == "OPH" & !is.na(phylum), .N, by = .(phylum)][order(-N)]$phylum[1:10]
c <- pd[primer == "UM" & !is.na(phylum), .N, by = .(phylum)][order(-N)]$phylum[1:10]
# get union of top phyla
phylist <- Reduce(union, list(a, b, c))
# reduce plotting data to most abundant phyla
pd <- pd[phylum %in% phylist]

# set colors for barplot
pal <- "Paired"
max <- brewer.pal.info[which(rownames(brewer.pal.info) == pal),]$maxcolors
cols <- colorRampPalette(brewer.pal(max, pal))(length(phylist))
# stacked barpchart
ggplot(pd, aes(fill = phylum, y = N, x = site)) + 
  geom_bar(position= "stack", stat = "identity", width = 0.5) + 
  scale_fill_manual(values = cols) + 
  labs(fill = "Phylum", y = "Number of taxa", x = "Sample") + 
  ylim(c(0, 680)) +
  facet_grid(. ~ primer, scales = "free_x", space = "free") +
  theme_bw() + theme(panel.grid = element_blank())






# plot Manila Bay
# packages
library(sf)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(gridExtra)

# https://r-spatial.github.io/sf/articles/sf1.html
# https://rstudio-pubs-static.s3.amazonaws.com/544700_a4314bd8b2d04d66b67e624dc0ccca7d.html

# set working directory
setwd("D:/Documents/NSRI 2/Datasets/")

# check layers of shape file
st_layers(dsn = "phl_adm_psa_namria_20200529_shp")

# read shape file
shp <- read_sf(dsn = "phl_adm_psa_namria_20200529_shp",
               layer = "phl_admbnda_adm3_psa_namria_20200529",
               quiet = FALSE)
# turn off spherical heometry
# https://github.com/r-spatial/sf/issues/1902
sf_use_s2(FALSE)
# set limits of bounding box
bb <- data.frame(xmin = 120.65, ymin = 14.4, 
                 xmax = 121.03, ymax = 14.8)
# crop shape file
x <- st_crop(shp, 
             xmin = bb$xmin, ymin = bb$ymin, 
             xmax = bb$xmax, ymax = bb$ymax)

# read-in coordinates of sites
readLines("D:/Documents/NGS/input/site_coordinates.csv")
coords <- read.csv("D:/Documents/NGS/input/site_coordinates.csv",
                   header = TRUE)[,3:5]
# get site # only
coords$site <- gsub("Site ", "", coords$site)

# plot Manila Bay map
map <- ggplot(x[, 3]) + 
  geom_point(data = coords, 
             aes(x = longitude, y = latitude),
             col = "red",
             size = 3) + 
  geom_text(data = coords,
            aes(x = longitude, y = latitude, label = site),
            size = 3, hjust = 1.5, vjust = -0.5) +
  geom_sf() +
  
  scale_x_continuous(breaks = seq(120.6, 121, by = 0.1)) + 
  scale_y_continuous(position = "right") +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text.y.right = element_text(color = "red"),
        axis.text.y.left = element_blank())

# stacked barpchart to add
################################
# read-in data
df <- readRDS("D:/Documents/NGS/out/Manila_Bay/transients/compiled_kraken_stdout_per_read.RDS")

# cleaning of dataset if blanks and extraction negatives are to be excluded
# subset df (remove unclassified reads) and restrict to COI samples since only those have blanks and negatives
sub <- df[taxid != "0" & primer == "UM"]

# get number of reads per taxon per replicate (samples)
samp <- sub[type == "sample", .(reads = .N), by = .(site, primer, taxid, taxon, taxon_rank_simp, superkingdom, kingdom, phylum, class, order, family, genus)]
# add corresponding blank per sample
samp$blank <- ifelse(grepl("^1|^2", samp$site), "b1B", 
                     ifelse(grepl("^8|^9", samp$site), "b8A", 
                            ifelse(grepl("^3|^4|^5", samp$site), "b5A", NA)))
# check assigned blanks
unique(samp[,c("site", "blank")])

# get number of reads per taxon per replicate (extraction negative)
neg <- sub[type == "negative", .(reads = .N), by = .(site, taxid)]
# get number of reads per taxon per replicate (blanks)
blanks <- sub[type == "blank", .(reads = .N), by = .(site, taxid)]

# get number of extraction negative reads that will be subtracted
samp$reads_neg <- neg$reads[match(samp$taxid, neg$taxid)]
# which taxa will be removed (those with reads less than or equal to zero)
samp[reads - reads_neg <= 0]

# list of replicates
reps <- unique(samp$site)
# add empty column
samp$reads_blank <- 0
# fill-in number of reads (from blanks) that will be subtracted
for(i in 1:length(reps)) {
  # subset per replicate
  s <- samp[site == reps[i]]
  # subset blanks to corresponding blank of replicate
  b <- blanks[site == unique(s$blank)]
  
  # match number of blank reads
  s$reads_blank <- b$reads[match(s$taxid, b$taxid)]
  samp[samp$site == reps[i]]$reads_blank <- s$reads_blank
}

# replace all NAs with 0
samp$reads_neg[is.na(samp$reads_neg)] <- 0
samp$reads_blank[is.na(samp$reads_blank)] <- 0
# get final reads excluding blanks and negative
samp$reads_fin <- samp$reads - samp$reads_neg - samp$reads_blank

# check number of taxids that will be removed
a <- samp[reads_blank > 0 | reads_neg > 0, c("site", "taxid", "reads", "reads_neg", "reads_blank", "reads_fin")]
a 
table(a$site) # taxa per replicate
table(a[reads_fin > 0]$site) # taxa retained per replicate
table(a[reads_fin <= 0]$site) # taxa to be removed per replicate

# format UM reads to be retained for merging with reads from other primers
temp <- samp[reads_fin > 0, -(13:16)]
# UM reads to be merged with other primers if blanks and negatives will not be removed:
# temp <- sub[type == "sample", .(reads_fin = .N), by = .(site, primer, taxid, taxon, taxon_rank_simp, superkingdom, kingdom, phylum, class, order, family, genus)]

# get reads from other primers
sub <- df[taxid != "0" & primer != "UM"]
# get number of reads per taxon per replicate (samples)
oth <- sub[type == "sample", .(reads_fin = .N), by = .(site, primer, taxid, taxon, taxon_rank_simp, superkingdom, kingdom, phylum, class, order, family, genus)]

# bind reads from UM with untouched reads from 16S and OPH
data <- rbind(oth, temp)
# consider only taxa identified until the class level
# data <- data[taxon_rank_simp == "class"]
# add site number
data$siteno <- gsub("[[:alpha:]]", "", data$site)

# get number of taxa per phylum (per siteno)
pd <- data[, .N, by = .(siteno, phylum, primer)]
# get the top most common phyla per primer
a <- pd[primer == "16S" & !is.na(phylum), .N, by = .(phylum)][order(-N)]$phylum[1:8]
b <- pd[primer == "OPH" & !is.na(phylum), .N, by = .(phylum)][order(-N)]$phylum[1:8]
c <- pd[primer == "UM" & !is.na(phylum), .N, by = .(phylum)][order(-N)]$phylum[1:8]
# get union of top phyla
phylist <- Reduce(union, list(a, b, c))
# reduce plotting data to most abundant phyla
pd <- pd[phylum %in% phylist]

# set order of sites
pd$siteno <- factor(pd$siteno, levels = c("6", "4", "7", "3", "5", "9", "8", "2", "1"))
# add location
pd$location <- ifelse(pd$siteno %in% c("6", "4", "7", "3", "5"), "north", "south")

# set colors for barplot
pal <- "Paired"
max <- brewer.pal.info[which(rownames(brewer.pal.info) == pal),]$maxcolors
cols <- colorRampPalette(brewer.pal(max, pal))(length(phylist))
# stacked barpchart
p <- ggplot(pd, aes(fill = phylum, y = N, x = primer)) + 
  geom_bar(position= "stack", stat = "identity", width = 0.5) + 
  scale_fill_manual(values = cols) + 
  labs(x = "", fill = "Phylum", y = "Number of taxa") + 
  ylim(c(0, 1090)) +
  facet_wrap( ~ siteno, nrow = 2, labeller = ) +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        legend.position = "left", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust = 1))
################################
map + inset_element(p, -.6, .1, .45, .9, 
                    clip = TRUE, 
                    align_to = "plot")
# grid.arrange(p, map, nrow = 1, ncol = 2)






# SCRATCH

# read-in data
df <- readRDS("D:/Documents/NGS/out/Manila_Bay/transients/compiled_kraken_stdout_per_read.RDS")
# check species in blanks
temp <- df[df$type == "blank" | df$type == "negative",]
temp <- df[df$type == "negative",]
# tabulation of rank
table(temp$taxon_rank_simp)
# species in blanks
unique(temp[taxon_rank_simp == "species"]$taxon)
sub <- unique(temp[taxon_rank_simp == "species", 10:21])
sub[phylum == "Echinodermata"]


# add site number to df
df$siteno <- gsub("\\D+", "", df$site)
unique(df[,c("site", "siteno")])
# subset df
temp <- df[, .(reads = .N), by = .(site, primer, taxid)]