#!/usr/bin/env python3

# setup ================================================================================
# import libraries
import numpy as np
import pyBigWig
import pandas as pd
import sys

# get command line arguments
in_fn = snakemake.input[0]
out_fn = snakemake.output[0]


# calculate genome-wide mean and sd ====================================================
# Get list of chromosomes to iterate over
bw = pyBigWig.open(in_fn)
bw_chroms = bw.chroms()

# calculate genome-wide mean
bw_header = bw.header()
total_bases = bw_header["nBasesCovered"]
gw_mean = bw_header["sumData"] / total_bases

# calculate genome wide SD
ssd = 0
print("calculating genome-wide standard deviation")
for chrom, n_bases in bw_chroms.items():
    #print(f"calculating sum of sqaured deviations for chromosome: {chrom}")
    chr_values = bw.values(chrom, 0, n_bases, numpy=True)
    chr_values = chr_values[~np.isnan(chr_values)]
    ssd += np.sum((chr_values - gw_mean)**2)

gw_sd = np.sqrt((ssd / (total_bases - 1)))

print("genome-wide mean",gw_mean)
print("genome-wide standard deviation",gw_sd)

# write z-score normalized values to output file ========================================
# open output file and add header
out_bw = pyBigWig.open(out_fn, "w")
out_bw.addHeader(list(bw_chroms.items()))

# perform z-score normalization for each chromosome and write values to file
print("performing z-score normalization on bigwig file")
for chrom, n_bases in bw_chroms.items():

	# calculate z-score normalized values
	print(f"calculating z-score for chromosome: {chrom}")
	chr_data = pd.DataFrame(bw.intervals(chrom), columns =['start', 'end', 'score'])
	chr_data["zscore"] = (chr_data['score'] - gw_mean) / gw_sd
	chr_data['chromosome'] = chrom
    
	# write output values to file
	out_bw.addEntries(
        list(chr_data["chromosome"]), 
        list(chr_data["start"]),
        ends = list(chr_data["end"]),
        values = list(chr_data["zscore"])
        )
   

# close bw files
bw.close()
out_bw.close()