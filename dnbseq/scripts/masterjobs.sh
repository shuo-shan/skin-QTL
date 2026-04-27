#!/bin/bash

ls /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/data/bam/softlinks | tr '\t' '\n' | cut -d"_" -f2,3 > temp.id.txt
while read x;do
  sed "s/id=/id=$x/g" job.sh > temp.job.sh
  bsub -o /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/scripts/reports/esat.$x.output.txt -J $x < temp.job.sh
  rm temp.job.sh
done < temp.id.txt
rm temp.id.txt
