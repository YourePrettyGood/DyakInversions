#!/usr/bin/env Rscript
#Define prefix, Kmax, and BEDfile
options <- commandArgs(trailingOnly=TRUE);
prefix <- options[1];
Kmax <- options[2];
BEDfile <- options[3];

cat(paste0("Using input and output file prefix: ", prefix, "\n"))
cat(paste0("Using maximum K of: ", Kmax, "\n"))
cat(paste0("Using windowing BED file: ", BEDfile, "\n"))

#Read in the BED:
window_bed <- read.table(BEDfile, colClasses=c("character", "integer", "integer"), col.names=c("Scaffold", "Start", "End"), header=FALSE);

#Now create the list to store Qmatrices:
Q_hats <- list();
num_skipped <- 0;
for (h in seq(1, nrow(window_bed))) {
   arm <- window_bed[h, "Scaffold"];
   #Declare a list if the scaffold hasn't been seen before:
   if (! arm %in% names(Q_hats)) {
      Q_hats[[arm]] <- list();
   }
   #Move the window counter forward:
   window_id <- length(Q_hats[[arm]])+1;
   #Read in the Q matrix for this window:
   #Be sure to skip the window if its Q matrix doesn't exist -- indicates ALStructure failure:
   if (!file.exists(paste0("QMatrices/", prefix, "_window", h, "_Qmatrix.tsv"))) {
      cat(paste("Q matrix missing for window", h, "-- did ALStructure fail for it?\n"));
      num_skipped <- num_skipped + 1;
      next;
   }
   Qmatrix <- read.table(paste0("QMatrices/", prefix, "_window", h, "_Qmatrix.tsv"), header=TRUE, row.names=1, stringsAsFactors=TRUE)
   #Determine the number of ancestral demes estimated for this window:
   K <- ncol(Qmatrix)
   cat(paste("Window", h, "has", K, "ancestral demes\n"))
   #Add some extra columns to the Q matrix to help with plotting:
   Qdf <- data.frame(Window=rep(window_id, nrow(Qmatrix)), K=rep(K, nrow(Qmatrix)), Individual=rownames(Qmatrix));
   #Add in the deme columns to the data.frame:
   for (j in 1:Kmax) {
      if (j <= K) {
         #Fill out the deme columns, with slight relabeling of columns:
         Qdf[, paste0("Deme", j)] <- Qmatrix[, paste0("deme", j)];
      } else {
         #Zero out any deme columns beyond the current K:
         Qdf[, paste0("Deme", j)] <- 0.0;
      }
   }
   #Store the Q matrix for this window in the list:
   Q_hats[[arm]][[window_id]] <- Qdf;
}

#Indicate how many windows were skipped in total:
cat(paste("Skipped a total of", num_skipped, "windows\n"))

#Save the list into an Rdata file for retrieval by plotting code:
save(Q_hats, file=paste0(prefix, "_ALStructure_results.Rdata"));
