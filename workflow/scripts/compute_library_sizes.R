# setup ------------------------------------------------------------------------
library(tidyverse)
library(ShortRead)

# define input files -----------------------------------------------------------
fq_files <- snakemake@input
sample <- snakemake@wildcards[["sample"]]

# compute spike-in percentages and scaling factors -----------------------------
# get total read count for each sample
count <- 0
for (i in seq(fq_files)) {
  fq <- fq_files[i]
  message(paste("computing total reads for file:", sample))
  count <-  count + countFastq(fq)$records
}

# add total read counts to output table ----------------------------------------
library_sizes <- tibble(sample_name = sample, total_reads = count)

# write output table to file ---------------------------------------------------
library_sizes |> write_tsv(snakemake@output[[1]])