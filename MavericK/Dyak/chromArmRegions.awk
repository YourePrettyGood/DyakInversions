#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS="\t";
}
#This regex matches X, 2L, 2R, 3L, 3R, 4, and Y
# (and some more that may occur in other Drosophila, like
#  2, 3, XL, XR)
#Basically, we just want to subset out windows from a
# tab-separated position-based file for the major chromosome
# arms (e.g. BED, VCF, GFF, etc.)
$1~/^[234XY][LR]?$/{
   print;
}
