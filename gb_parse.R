# script for parsing GenBank flatfiles

# working directory: path of .gb files
setwd("D:/Documents/NGS/entrez/wormstaxlist_gb/compiled_gb")

# view possible formats
readLines("Aaptos_aaptos_ncbi.gb") # 1 COI sequence
readLines("Aaptos_papillata_ncbi.gb") # 2 COI sequences
readLines("Anemonia_viridis_ncbi.gb", 705) # complete genome

gb <- readLines("Aaptos_aaptos_ncbi.gb")
getseqInfo(gb)


getseqInfo <- function(genbankFile){
  temp <- grep(x = genbankFile, pattern = "ORIGIN|//")
  seqInfo <- genbankFile[(temp[1] + 1):(temp[2] - 1)]
  rm(temp)
  seqInfo <- paste(gsub(x = seqInfo, pattern = "[[:digit:]]| ",
                        replacement = ""),
                   collapse = "")
  return(seqInfo)
}

substr()