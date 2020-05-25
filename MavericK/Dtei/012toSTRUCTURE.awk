#!/bin/awk -f
BEGIN{
   OFS="\t";
}
{
   line1=$1;
   for (i=2;i<=NF;i++) {
      if ($i==0){
         line1=line1"\t0";
      }else if ($i==1){
         line1=line1"\t0";
      }else{
         line1=line1"\t1";
      };
   };
   print line1;
   line2=$1;
   for (i=2;i<=NF;i++) {
      if ($i==0){
         line2=line2"\t0";
      }else if ($i==1){
         line2=line2"\t1";
      }else{
         line2=line2"\t1";
      };
   };
   print line2;
}
