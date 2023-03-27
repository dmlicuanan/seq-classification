# load packages
library(xlsx)
library(data.table)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

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
# wide <- dcast(dt, sample ~ process, value.var = "reads")
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
# primer
dt$primer <- ifelse(grepl("UM", dt$sample), "UM", 
                    ifelse(grepl("16S", dt$sample), "16S", "OPH"))
# set (merged, unmerged, or no merging)
dt$set <- ifelse(grepl("^merged", dt$process), "Merged",
                 ifelse(grepl("^unmerged", dt$process), "Unmerged",
                        ifelse(grepl("^trimmed", dt$process), "No merging", "All")))

# remove reads after filtering (since classified sequences came straight from trimming, not filtering)
dt <- dt[!is.na(dt$stage), ]
# deciles of reads
quantile(dt$reads, probs=seq(0.1, 1, by = 0.1))

# plot for all, regardless of set
plot(x = dt$stage, dt$reads)

# limit y axis
plot(x = dt$stage, dt$reads, type = "n",
     xlab = "stage", ylab = "no. of reads")

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





# line graph for reads that remain after merging and trimming 
# proceed with dt data.table from above
# reads from filtering step are already excluded

# see unique values for variables
lapply(list(dt$process, dt$stage, dt$set), unique)
# check values of set
unique(dt[, c("set", "process")])

# remove reads which did not undergo merging
m <- dt[set != "No merging"]
# remove blanks and negatives
m <- m[!grepl("^[[:alpha:]]", m$sample), ]

# separate raw reads 
raw <- m[set == "All"]
# change set to "Unmerged"
raw$set <- "Unmerged"
# make "All" > "Merged"
m$set[grepl("^All", m$set)] <- "Merged"
# bind raw to m
m <- rbind(raw, m)

# add column for replicate name
m$rep <- sub("-.*", "", m$sample)
# add site number
m$site <- sub("(^[0-9]+).*$", "\\1", m$sample)

# make set a factor
m$set <- factor(m$set)

# subset of data to be plotted
pd <- m
# subset for main plot
pd <- m[set == "Merged"]
pd <- m[set == "Unmerged"] # add linetype = "dashed" outside geom_line's aes
# subset for inset plot
pd <- m[set == "Merged" & primer != "UM"]

# plot of number of reads remaining from samples at each stage of processing
ggplot(pd, aes(x = stage, 
               y = reads, 
               group = interaction(sample, set),
               linetype = set)) + 
  geom_point(aes(col = primer)) + 
  geom_line(aes(col = primer)) + 
  scale_colour_brewer(palette = "Set1") + 
  labs(x = "Process", 
       y = "Number of reads", 
       linetype = "", 
       col = "Primer") + 
  scale_x_continuous(labels = c("Sequencing", "Merging", 
                                "Trimming & filtering", 
                                "Classification")) +
  theme_minimal() +  
  theme(panel.grid.minor.x = element_blank()) 

# merged only (number of reads per stage)
# main plot
main <- ggplot(pd, aes(x = stage, 
              y = reads, 
              group = sample)) + 
  geom_point(aes(col = primer)) + 
  geom_line(aes(col = primer)) + 
  scale_colour_brewer(palette = "Set1") + 
  labs(x = "Process", 
       y = "Number of reads", 
       linetype = "", 
       col = "Primer") + 
  scale_x_continuous(labels = c("Sequencing", "Merging", 
                                "Trimming & filtering", 
                                "Classification"))  + 
  scale_y_continuous(labels = comma) +
  theme_classic() +  
  theme(panel.grid.minor.x = element_blank()) 

# inset plot
inset <- ggplot(pd, aes(x = stage, 
               y = reads, 
               group = sample)) + 
  geom_point(aes(col = primer)) + 
  geom_line(aes(col = primer)) + 
  scale_colour_brewer(palette = "Set1") + 
  scale_y_continuous(labels = comma) +
  labs(x = "Process", 
       y = "Number of reads", 
       linetype = "", 
       col = "Primer") +
  theme_classic() +  
  theme(panel.grid.minor.x = element_blank(), 
        axis.title = element_blank(), 
        axis.text.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x = element_blank(),
        legend.position = "none") 
inset
# combine main and inset plot
main + inset_element(inset, .6, .6, 1, 1, align_to = "plot")

# summary of number of reads remaining
pd <- m[,.(mean = mean(reads), sd = sd(reads)), by = .(stage, primer, set)]
pd[order(set, primer, stage)]
# plot means
ggplot(pd, aes(x = stage, 
               y = mean, 
               linetype = set)) + 
  geom_point(aes(col = primer)) + 
  geom_line(aes(col = primer)) + 
  scale_colour_brewer(palette = "Set1") + 
  labs(x = "Process", 
       y = "Number of reads", 
       linetype = "", 
       col = "Primer",
       title = "Mean number of reads") + 
  scale_x_continuous(labels = c("Sequencing", "Merging", 
                                "Trimming & filtering", 
                                "Classification")) +
  theme_minimal() +  
  theme(panel.grid.minor.x = element_blank())
