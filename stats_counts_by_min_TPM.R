#!/usr/bin/env Rscript


# Title: stats_counts_by_min_TPM.R
# Date: 2018-01-11
# Author: Adam Nunn
#
# Description:
#    This R script carries out a linear regression analysis on the first
#    100 values for number of expressed features corresponding to
#    incremental increases in the TPM threshold, and plots a graph. The
#    script takes a single input file: 1) the table of values generated
#    by the trinity script "count_matrix_features_given_MIN_TPM_threshold.pl"
#    (INPUT_TPM_FILE). The script will output the graph with a user-specified
#    upper limit on the y-axis (UPPER_LIMIT) and produce the file
#    "matrix.TPM.not_cross_norm.counts_by_min_TPM.pdf" in the current
#    directory. The script will also print the y-intercept to STDOUT.
#
# Usage:
#     Rscript stats_counts_by_min_TPM.R UPPER_LIMIT INPUT_TPM_FILE
# eg. Rscript stats_counts_by_min_TPM.R 25000 matrix.TPM.not_cross_norm.counts_by_min_TPM
#

args <- commandArgs(trailingOnly=T)
data <- read.table(args[2], header=T)

# extract the data between 10 TPM and 100 TPM
filt_data = data[data[,1] > -100 & data[,1] < -10,] 
# perform a linear regression on this filtered subset of the data
fit = lm(filt_data[,2] ~ filt_data[,1])
print(fit)
 
# add the linear regression line to the plot
pdf("matrix.TPM.not_cross_norm.counts_by_min_TPM.pdf")
plot(data, xlim=c(-100,0), ylim=c(0,as.integer(args[1])), t='b')
abline(fit, col='green', lwd=3)
dev.off()

