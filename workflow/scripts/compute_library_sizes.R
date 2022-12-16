# setup ------------------------------------------------------------------------
library(tidyverse)
library(ShortRead)

# define input files -----------------------------------------------------------
fq_files <- snakemake@input[[1]]
sample <- snakemake@wildcards[["sample"]]

# compute spike-in percentages and scaling factors -----------------------------
# get total read count for each sample
  fq <- fq_files[1]
  message(paste("computing total reads for file:", sample))
  count <-  countFastq(fq)$records

# add total read counts to output table ----------------------------------------
library_sizes <- tibble(sample_name = sample, total_reads = count)

# write output table to file ---------------------------------------------------
library_sizes |> write_tsv(snakemake@output[[1]])