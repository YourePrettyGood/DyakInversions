#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
#If we don't pass in genomewide as a non-zero value, output per-scaffold:
   if (length(genomewide)==0) {
      genomewide=0;
   };
#Account for weights in column 4 if indicated, otherwise consider it omit:
   if (length(weighted)==0) {
      weighted=0;
   };
#Don't use weights for weighted average, but filter on this infimum if set:
   if (length(infimum)==0) {
      infimum=-1;
   }
#Use a different column than the default 3 for stat:
   if (length(statcol)==0) {
      statcol=3;
   }
}
{
#Only include sites without the omission column if that column is present:
   if (weighted == 0 && (NF == 3 || $4 == 0)) {
      if (statcol > 3) {
         print "Nonsensical statcol "statcol" supplied. Bugging out." > "/dev/stderr";
         exit 1;
      }
      stat[$1]+=$statcol;
      count[$1]+=1;
   } else if (weighted == 1) {
#If weighted is selected, but we don't have a weight column, revert:
      if (weighted != 0 && NF > 3) {
         if (statcol > NF) {
            print "Nonsensical statcol "statcol" supplied. Bugging out." > "/dev/stderr";
            exit 2;
         }
         if (infimum != -1) {
            if ($4 > infimum) {
               stat[$1]+=$statcol;
               count[$1]+=1;
            }
#Use the weighted average if infimum isn't set:
         } else {
            stat[$1]+=$statcol*$4;
            count[$1]+=$4;
         }
      } else {
         if (statcol > NF) {
            print "Nonsensical statcol "statcol" supplied. Bugging out." > "/dev/stderr";
            exit 3;
         }
         stat[$1]+=$statcol;
         count[$1]+=1;
      };
   };
}
END{
   if (genomewide != 0) {
      for (scaf in stat) {
         stat_gw+=stat[scaf];
         count_gw+=count[scaf];
      };
      if (count_gw > 0) {
         print "Genome-wide", stat_gw/count_gw;
      } else {
         print "Genome-wide", "NA";
      };
   } else {
      for (scaf in stat) {
         if (count[scaf] > 0) {
            print scaf,stat[scaf]/count[scaf];
         } else {
            print scaf,"NA";
         };
      };
   };
}
