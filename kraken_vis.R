# location of kraken files
setwd("D:/Documents/NGS/out/Manila_Bay/kraken/")

# load packages
library(data.table)
library(stringr)

# preview kraken standard output and report files
# kraken standard output and report files
f1 <- list.files(pattern = ".kraken$")
f2 <- list.files(pattern = ".krakenReport$")
# index of file to be read
i <- 1

# preview file
readLines(f1[i], 10)
readLines(f2[i], 10)

# load file
d1 <- fread(f1[i], sep = "\t", col.names = c("code", "seqid", "taxon", "rdlen", "lcamap"))
d2 <- fread(f2[i], sep = "\t", col.names = c("percroot", "ctroot", "ctdirect", "rank", "taxid", "sciname"))

# nodes.dmp
nodes <- fread("D:/Documents/NGS/entrez/taxonomy/nodes.dmp", sep = "|", drop = c(7:10,12:14), col.names = c("taxid", "parent", "rank", "emblcode", "divid", "divflag", "hidflag"))
nodes <- as.data.table(lapply(nodes, trimws))




# list of files to be compiled
# base pattern on primer, sites/blanks/negatives
files <- list.files(pattern = "^[1-9].*_merged_trimmed.kraken$")

# empty list
df <- vector(mode = "list", length = length(files))
# read each standard output file
for(i in 1:length(files)) {
  temp <- fread(files[i], sep = "\t", col.names = c("code", "seqid", "taxon", "rdlen", "lcamap"))
  temp <- temp[, .N, by = .(taxon, code)]
  temp$site <- sub("-.*", "", files[i])
  temp$taxid <- gsub(".*\\(taxid (.+)\\)", "\\1", temp$taxon)
  df[[i]] <- temp
}

# collapse all stdout in one list
df <- rbindlist(df)
# make columns taxa, and rows the sites
mat <- dcast(df, taxid ~ site, value.var = "N")
# replace all NAs with 0
mat[is.na(mat)] <- 0
# make 1st column row names
rownames(mat) <- mat$taxid
 <- samp[,-1]


# https://stackoverflow.com/questions/9617348/reshape-three-column-data-frame-to-matrix-long-to-wide-format

