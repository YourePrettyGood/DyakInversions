#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   #If we don't pass in infimum as in [0, 1), set it to 0:
   if (length(infimum)==0 || infimum < 0 || infimum >= 1.0) {
      infimum=0;
   };
   #If no column number is passed, use the third column:
   if (length(statcol)==0) {
      statcol=3;
   }
   #Use this flag to trigger segregating sites count
   # i.e. sum up if statcol value is > 0, not the actual stat
   # Given pi or H in the statcol, the sum is then S
   if (length(segsites)==0) {
      segsites=0;
   }
   n_intervals=0;
}
#First file is a BED of intervals to evaluate:
FNR==NR{
   n_intervals+=1;
   intervals["scaf",n_intervals]=$1;
   intervals["start",n_intervals]=$2+1;
   intervals["end",n_intervals]=$3;
   stat[n_intervals]=0;
   count[n_intervals]=0;
}
#Second file is the statistics file
#We will be inefficient and search through all intervals each time
FNR<NR&&FNR==1{
   print "Header was: " $0 > "/dev/stderr";
}
FNR<NR&&FNR>1{
   #Error out if there aren't enough columns:
   if (NF < 4) {
      print "Not enough columns, we expect at least 4 and there were "NF". Bugging out." > "/dev/stderr";
      exit 1;
   }
   #Iterate through the intervals:
   for (i=1; i<=n_intervals; i++) {
      if ($1 == intervals["scaf",i] && $2 >= intervals["start",i] && $2 <= intervals["end",i]) {
         if ($4 > infimum) {
            if (segsites) {
               stat[i]+=$statcol > 0 ? 1 : 0;
            } else {
               stat[i]+=$statcol;
            }
            count[i]+=1;
         }
      }
   }
}
END{
   for (i=1; i<=n_intervals; i++) {
      print intervals["scaf",i], intervals["start",i]-1, intervals["end",i], stat[i], count[i];
   }
}
