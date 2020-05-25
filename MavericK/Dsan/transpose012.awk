#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   #Keep track of how many sites so final iteration is known order:
   n_sites=0;
}
#Grab the individual IDs, ignore SNP column:
NR==1{
   for (i=2; i<=NF; i++) {
      indivs[i-1]=$i;
   }
   n_indivs=NF-1;
}
#Now read in all the genotypes and store
NR>1{
   n_sites+=1;
   if (NF != n_indivs + 1) {
      print "Somehow the matrix isn't rectangular at line "NR", bugging out." > "/dev/stderr";
      exit 1;
   }
   for (i=2; i<=NF; i++) {
      genotypes[indivs[i-1],n_sites]=$i;
   }
}
END{
   for (i=1; i<=n_indivs; i++) {
      output=indivs[i];
      for (j=1; j<=n_sites; j++) {
         output=output"\t"genotypes[indivs[i],j];
      }
      print output;
   }
}
