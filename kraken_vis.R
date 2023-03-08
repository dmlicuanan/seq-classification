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
# restrict entries to those with replicates (sample sequenced using uniminibarcode primer)
# remove site 4 since it has no replicate
# remove unclassified reads
sub <- df[primer == "UM" & type == "sample" & site != "4A" & code == "C"]
# summarize number of reads by site, etc (so that there is one row for each taxon in a replicate)
sub <- sub[, .(reads = .N), by = .(site, taxid, taxon, taxon_rank_simp, superkingdom, kingdom, phylum, class, order, family, genus)]
# add site number
sub$siteno <- gsub("\\D+", "", sub$site)
unique(sub[,c("site", "siteno")])

# create empty list
names <- c("site", "total", "mean", "sem", "common", "shared", "unique")
l <- sapply(names, function(x) NULL)

# add names of sites with replicates
l$site <- unique(sub$siteno)

for(i in 1:length(l$site)) {
  # create subset for site
  temp <- sub[siteno == l$site[i]]
  # set taxonomic rank to be used for making counts 
  temp$x <- temp$taxid
  
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
  }
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