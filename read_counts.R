# load packages
library(xlsx)
library(data.table)
library(ggplot2)

# read-in raw read counts
d <- xlsx::read.xlsx(file = "D:/Documents/NGS/out/Manila_Bay/transients/line_counts.xlsx",
                     sheetName = "# lines", colIndex = 1:13)

# convert all columns to numeric
d[,2:13] <- lapply(d[,2:13], as.numeric)

# except for columns with "classified", where 1 line is equivalent to 1 read,
# numeric columns need to be divided by 4 
# (since counts here are based on fastq files, where 1 read is equivalent to 4 lines)
# indices of the columns whose values need to be divided by 4:
ind <- grep(pattern = "classified|sample", x = colnames(d), invert = TRUE)
# get actual read counts from line counts
d[,ind] <- d[,ind]/4

# convert to data.table
dt <- data.table(d)
# convert from wide to long
dt <- melt(dt, id.vars = "sample", measure.vars = colnames(d)[-1])

# rename variables
colnames(dt) <- c("sample", "process", "reads")

# create read_counts.xlsx
# convert from long to wide
wide <- dcast(dt, sample ~ process, value.var = "reads")
# write.xlsx(wide, file = "D:/Documents/NGS/out/Manila_Bay/transients/read_counts.xlsx", row.names = FALSE )


# add additional columns
# stage (excluding filtering stage since most samples are left unfiltered)
dt$stage <- ifelse(dt$process == "raw", 1, 
                   ifelse(grepl("merged$", dt$process), 2, 
                          ifelse(grepl("trimmed$", dt$process), 3,
                                 ifelse(grepl("classified$", dt$process), 4, NA))))
# 1. Number of raw reads
# 2. Number of reads after merging
# 3. Number of reads after trimming
# 4. Number of reads after classification

# gene (based on primer)
dt$gene <- ifelse(grepl("UM", dt$sample), "COI", "16S")
# set (merged, unmerged, or no merging)
dt$set <- ifelse(grepl("^merged", dt$process), "Merged",
                 ifelse(grepl("^unmerged", dt$process), "Unmerged",
                        ifelse(grepl("^trimmed", dt$process), "No merging", "All")))

# remove entries where reads = NA (those with no result after filtering)
dt <- dt[!is.na(dt$reads), ]
# deciles of reads
quantile(dt$reads, probs=seq(0.1, 1, by = 0.1))

# plot for all, regardless of set
plot(x = dt$stage, dt$reads)

# limit y axis
plot(x = dt$stage, dt$reads, type = "n",
     xlab = "stage", ylab = "no. of reads", ylim = c(0,500000))

# samples
samp <- unique(dt$sample)

# plot for merged
# subset based on Merged
sub <- dt[(dt$set == "All" | dt$set == "Merged") & !is.na(dt$stage),]
# add lines per sample
for(i in 1:34) {
  temp <- sub[sub$sample == samp[i],]
  points(x = temp$stage, y = temp$reads, type = "o", col = "green")
}

# plot for reads that did not undergo merging
# subset based on No merging
sub <- dt[(dt$set == "All" | dt$set == "No merging") & !is.na(dt$stage),]
# add lines per sample
for(i in 1:34) {
  temp <- sub[sub$sample == samp[i],]
  points(x = temp$stage, y = temp$reads, type = "o", col = "blue")
}

# legend
legend(x = "topright", lty = 1, col= c("green","blue"),
       legend=c("Merged", "Trimmed only"), 
       bty = "n")

# proceed with using reads that undergo Merging
# see read_counts.xlsx 