#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
}
#Set up the header line:
NR==1{
   #Replace the [1]CHROM:POS column header with "SNP", since it's a SNP ID:
   $1="SNP";
   for (i=2; i<=NF; i++) {
      #Remove the [n] from each column header, and eliminate the colon:
      split($i, colid, "[]:]");
      $i=colid[2];
   };
   print $0;
}
NR>1{
   for (i=2; i<=NF; i++) {
      #Split the genotype into alleles (allow for unphased or phased genotype):
      split($i, gt, "[/|]");
      #Add the allele IDs (diploid-only) to get the number of alt alleles:
      $i=gt[1]+gt[2];
   };
   print $0;
}
