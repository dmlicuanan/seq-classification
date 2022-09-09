# script for formatting bold .fastas to kraken format

# readFasta()
# Reads .fas file and returns data frame
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

# writeFasta()
# Write .fas file from a vector of IDs (accession numbers info etc.) and sequences
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

