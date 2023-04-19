# !/usr/bin/env Rscript

library(deseq2)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

# Load the raw counts data into a DGEList object
#counts <- read.table(infile, header=TRUE, row.names=1, sep='\t')
#counts <- as.matrix(counts[, -1]) # remove first column of gene names
