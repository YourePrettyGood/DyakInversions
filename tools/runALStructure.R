#!/usr/bin/env Rscript
#Arguments:
#1) Path to input TSV, which consists of rows as SNPs, columns as individuals,
#    and there is an extra row and column with row and column names, respectively
#2) Path to desired output TSV, which is the labeled transpose of the Q_hat matrix
#    That is, columns are ancestral demes, rows are individuals.
#    We store the transpose because this matches the format of the Q matrix output
#    by MavericK.

options <- commandArgs(trailingOnly=TRUE)

#Load the necessary libraries (expects them to already be installed):
library(tidyverse)
library(alstructure)

#Read the input parameters from the command line:
input_012_tsv <- options[1]
output_Qhat_tsv <- options[2]
output_Phat_tsv <- options[3]

#Read in the 012 genotype data:
filtered_snps <- as.matrix(read.table(input_012_tsv, header=TRUE, row.names=1, stringsAsFactors=FALSE))

#Print the dimensions of the input:
cat(paste(ncol(filtered_snps), "individuals at", nrow(filtered_snps), "SNPs for", input_012_tsv))

#Run ALStructure (defaults for now):
filtered_ALStructure <- tryCatch(alstructure(filtered_snps),
                                 error=function(x) {"unsolvable"})
#Check if ALStructure failed, and skip output if it did:
if (class(filtered_ALStructure) != "character") {
   #Extract Q_hat and label the individuals and demes, then transpose
   # to match format of Q matrices from MavericK:
   filtered_Q_hat <- filtered_ALStructure$Q_hat
   colnames(filtered_Q_hat) <- colnames(filtered_snps)
   rownames(filtered_Q_hat) <- paste0("deme", 1:nrow(filtered_Q_hat))
   #We can prettify the names later to remove _[species] or [species]_
   # but we leave this code generic by skipping the prettification
   Q_hat_T <- t(filtered_Q_hat)

   #Print out the number of ancestral demes estimated:
   cat(paste(ncol(Q_hat_T), "ancestral demes for window", input_012_tsv))

   #Write out Q_hat:
   write.table(Q_hat_T, file=output_Qhat_tsv, quote=FALSE, sep="\t")

   #Briefly output labeled P_hat:
   filtered_P_hat <- as.data.frame(filtered_ALStructure$P_hat)
   colnames(filtered_P_hat) <- paste0("deme", 1:ncol(filtered_P_hat))
   filtered_P_hat$scaf <- sapply(rownames(filtered_P_hat), function(x) {strsplit(x, ":", fixed=TRUE)[[1]][1]})
   filtered_P_hat$pos <- sapply(rownames(filtered_P_hat), function(x) {strsplit(x, ":", fixed=TRUE)[[1]][2]})

   write.table(filtered_P_hat, file=output_Phat_tsv, quote=FALSE, sep="\t")
   
   #Clean up:
   rm(filtered_Q_hat, Q_hat_T, filtered_P_hat)
} else {
   cat(paste("Skipping output of Q matrix for window", input_012_tsv, "due to ALStructure failure"))
}
#Clean up:
rm(filtered_snps, filtered_ALStructure)