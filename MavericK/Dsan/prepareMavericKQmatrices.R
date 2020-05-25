#!/usr/bin/env Rscript
#Define prefix, Kmax, BEDfile, and sufwidth
options <- commandArgs(trailingOnly=TRUE);
#Prefix for input and output files:
#e.g. input Qmatrix file would be paste0("QMatrices/", prefix, "_window", window_suffix, "_Qmatrix_perInd_K", K, ".csv")
prefix <- options[1];
#Maximum K assayed by MavericK:
Kmax <- options[2];
#BED file of windows on main chromosome arms:
BEDfile <- options[3];

cat(paste0("Using input and output file prefix: ", prefix, "\n"))
cat(paste0("Using maximum K of: ", Kmax, "\n"))
cat(paste0("Using windowing BED file: ", BEDfile, "\n"))

#Read in the BED:
window_bed <- read.table(BEDfile, colClasses=c("character", "integer", "integer"), col.names=c("Scaffold", "Start", "End"), header=FALSE);

#Determine the set of individuals for any missing windows:
for (h in seq(1, nrow(window_bed))) {
   if (file.exists(paste0("QMatrices/", prefix, "_window", h, "_Qmatrix_perInd_K1.csv"))) {
      QmatrixK1 <- read.csv(paste0("QMatrices/", prefix, "_window", h, "_Qmatrix_perInd_K1.csv"));
      break;
   }
}

#Now create the list to store Qmatrices:
MavericK <- list();
for (h in seq(1, nrow(window_bed))) {
   arm <- window_bed[h, "Scaffold"];
   #Declare a list if the scaffold hasn't been seen before:
   if (! arm %in% names(MavericK)) {
      MavericK[[arm]] <- list();
   }
   #Move the window counter forward:
   window_id <- length(MavericK[[arm]])+1;
   if (file.exists(paste0("Evidence/", prefix, "_window", h, "_Evidence.csv"))) {
      #Read in the log-likelihoods for each K:
      Kchoice <- read.csv(paste0("Evidence/", prefix, "_window", h, "_Evidence.csv"));
      #Select the K with the largest log-likelihood according to Thermodynamic Integration:
      K <- Kchoice[Kchoice$logEvidence_TI == max(Kchoice$logEvidence_TI), "K"];
      if (length(K) > 1) {
         cat(paste0("Found multiple best K for window ", h, "\n"))
         #Choose the first K, arbitrarily, since it should be the smallest K
         K <- K[1]
      }
      #Read in the Q matrix for the best K:
      Qmatrix <- read.csv(paste0("QMatrices/", prefix, "_window", h, "_Qmatrix_perInd_K", K, ".csv"));
      #Add some extra columns to the Q matrix to help with plotting:
      Qdf <- data.frame(Window=rep(window_id, nrow(Qmatrix)), K=rep(K, nrow(Qmatrix)), Individual=Qmatrix$label);
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
      MavericK[[arm]][[window_id]] <- Qdf;
   } else {
      Qdf <- data.frame(Window=rep(window_id, nrow(QmatrixK1)), K=rep(1, nrow(QmatrixK1)), Individual=QmatrixK1$label);
      for (j in 1:Kmax) {
         Qdf[, paste0("Deme", j)] <- 0.0;
      }
      MavericK[[arm]][[window_id]] <- Qdf;
      cat(paste0(arm, " window ", window_id, " aka window ", h, " not found, skipping.\n"))
   }
}

#Save the list into an Rdata file for retrieval by plotting code:
save(MavericK, file=paste0(prefix, "_MavericK_results.Rdata"));
