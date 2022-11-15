# reformat kraken output

# load packages
library(data.table)





# location of kraken files
setwd("D:/Documents/NGS/out/Manila_Bay/kraken/")

# preview kraken standard output and report files
# kraken standard output and report files
f1 <- list.files(pattern = ".kraken$")
f2 <- list.files(pattern = ".krakenReport$")
# index of file to be read
i <- 1

# preview file
readLines(f1[i], 10)
readLines(f2[i], 10)

# load file to see format
d1 <- fread(f1[i], sep = "\t", col.names = c("code", "seqid", "taxon", "rdlen", "lcamap"))
d2 <- fread(f2[i], sep = "\t", col.names = c("percroot", "ctroot", "ctdirect", "rank", "taxid", "sciname"))





# reformat nodes.dmp to table form (where species can be linked to all taxonomic information from species to superkingdom)
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
# first element are taxids in nodes.dmp in species level
x <- list(nodes[nodes$rank == "species",]$taxid) 
# add species as names of vector
# names(x[[1]]) <- rep("species", length(x[[1]]))

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
# convert from wide to long format, where 1st column is species taxid
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
dt$species <- names$name[match(dt$taxid, names$taxid)]
dt$parent_taxon <- names$name[match(dt$parent_taxid, names$taxid)]

# rearrange columns
dt <- dt[ , c(1,4,2,5,3)]
# save as RDS
# saveRDS(dt, "D:/Documents/NGS/entrez/taxonomy/taxonomy_long.RDS")
# save as .txt table
# write.table(dt, file = "D:/Documents/NGS/entrez/taxonomy/taxonomy_long.txt", sep = "\t", row.names = FALSE)







########### not sure if table should still be made wide given that some ranks do not hold unique values (e.g. "clade", "no rank")
# read RDS
dt <- readRDS("D:/Documents/NGS/entrez/taxonomy/taxonomy_long.RDS")

dcast(y, taxid + species ~ parent_rank, value.var = "parent_taxon")

# restrict table values to those unique ranks
x <- dt[dt$parent_rank != "no rank" & dt$parent_rank != "clade", ]
# wide format
wide <- dcast(x, taxid + species ~ parent_rank, value.var = c("parent_taxon", "parent_taxid"))
# arrange by taxid
wide <- wide[order(as.integer(taxid))]
# save .txt table
# write.table(wide, file = "D:/Documents/NGS/entrez/taxonomy/taxonomy_wide.txt", sep = "\t", row.names = FALSE, quote = TRUE)
# save RDS
# saveRDS(wide, "D:/Documents/NGS/entrez/taxonomy/taxonomy_wide.RDS")
###################################




# location of kraken files
setwd("D:/Documents/NGS/out/Manila_Bay/kraken/")
# list of kraken standard output to be compiled
# base pattern on primer, sites/blanks/negatives
files <- list.files(pattern = "^[1-9].*_merged_trimmed.kraken$")

# read-in files for adding other taxonomic information
# long taxonomy table
long <- readRDS("D:/Documents/NGS/entrez/taxonomy/taxonomy_long.RDS")

# create empty list
df <- vector(mode = "list", length = length(files))
# read each standard output file to fill-in list
for(i in 1:length(files)) {
  temp <- fread(files[i], sep = "\t", col.names = c("code", "seqid", "taxon", "rdlen", "lcamap"))
  # summarize counts per taxon
  temp <- temp[, .N, by = .(taxon, code)]
  # tag taxon with site name
  temp$site <- sub("-.*", "", files[i])
  # extract taxid
  temp$taxid <- gsub(".*\\(taxid (.+)\\)", "\\1", temp$taxon)
  df[[i]] <- temp
}

# collapse all kraken stdout in one list
df <- rbindlist(df)

# add rank 
df$rank <- long$parent_rank[match(df$taxid, long$parent_taxid)]

# add taxa information for kraken classification
taxget <- function(ktaxid, rank) {
  linker <- long[parent_taxid == ktaxid, ]$taxid[1]
  taxon <- long[taxid == linker & parent_rank == rank, ]$parent_taxon
  if (length(taxon) == 1) { return(taxon) } else { return(NA) }
}

# get list of unique taxids that taxonomic info are needed for
list <- data.table(taxid = unique(df$taxid))

for(i in 1:length(list$taxid)) {
  list$phylum[i] <- taxget(list$taxid[i], "phylum")
  list$class[i] <- taxget(list$taxid[i], "class")
  list$order[i] <- taxget(list$taxid[i], "order")
  list$family[i] <- taxget(list$taxid[i], "family")
  print(i)
}

# append to df table
df$phylum <- list$phylum[match(df$taxid, list$taxid)]
df$class <- list$class[match(df$taxid, list$taxid)]
df$order <- list$order[match(df$taxid, list$taxid)]
df$family <- list$family[match(df$taxid, list$taxid)]











# make columns taxa, and rows the sites
dt <- dcast(df, site ~ taxid, value.var = "N")
# replace all NAs with 0
dt[is.na(dt)] <- 0

# create matrix  
mat <- as.matrix(dt[,-1], rownames = dt$site)
# replace all non-zeros with 1
mat[mat > 0] <- 1

d <- vegdist(mat)

# hclust
h <- hclust(d)
plot(h)

m <- metaMDS(d, display = "site", )
plot(m, type = "t")
ordihull(m, rownames(mat))
