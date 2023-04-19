# !/usr/bin/env Rscript

# Description: Normalise the raw counts using the TMM method
# Author: Nienke van Unen
# Date: 14-04-2023

# Usage: Rscript normalise_counts.R <infile> <outfile>
# Arguments:
#  infile: path to the raw counts file
#  outfile: path to the output file
# Output: a file with the normalised counts
# Dependencies: edgeR
# Notes: This script is called by the Snakefile
# Example: Rscript normalise_counts.R raw_counts.txt normalised_counts.txt
# Version: 1.0
# History: 14-04-2023: Initial version

library(edgeR)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

# Load the raw counts data into a DGEList object
counts <- read.table(infile, header=TRUE, row.names=1, sep='\t')
#counts <- as.matrix(counts[, -1]) # remove first column of gene names

dge <- DGEList(counts)

# Perform TMM normalization
dge <- calcNormFactors(dge, method='TMM')

# Extract the normalized counts
norm_counts <- cpm(dge)

# Write the normalized counts to a file
write.table(norm_counts, file=outfile, sep='\t', quote=FALSE)