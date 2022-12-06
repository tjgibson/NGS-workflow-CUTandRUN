# setup ------------------------------------------------------------------------
library(tidyverse)
library(ShortRead)



# read in nucleosome barcodes --------------------------------------------------
nuc_barcodes <- read_csv("config/epicypher_met_panel_nuc_barcodes.csv")

# get input fastq files --------------------------------------------------------
input_files <- list.files("data/raw_reads/", full.names = TRUE)
col_names <- basename(input_files) |> str_replace("_001.fastq.gz", "")

nuc_barcodes[,col_names] <- 0

# read in raw read sequences ---------------------------------------------------
chunk_size <- 1e6

for (i in seq(col_names)) {
  fq <- input_files[i]
  
  strm <- FastqStreamer(fq, readerBlockSize = chunk_size)
  on.exit(close(stream))
  
  message(paste("starting barcode search for file:", fq))  
  
  record_count <- 0
  repeat {
    record_count <- record_count + chunk_size
    message(paste("reads processed:", record_count))
    chunk <- yield(strm)
    if (length(chunk) == 0)
      break
    ## process chunk
    for (j in seq(nrow(nuc_barcodes))) {
      chunk_count <- vcountPattern(pull(nuc_barcodes[j,"sequence"]), sread(chunk)) |> sum()
      nuc_barcodes[j,col_names[i]] <- nuc_barcodes[j,col_names[i]] + chunk_count
      
    }
  }
  print(nuc_barcodes)
}

nuc_barcodes |> write_tsv("results/spike_in_nuc_counts_test.txt")
