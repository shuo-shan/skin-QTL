#!/bin/bash
### script for CELseq2 data processing
# change read name to have :barcode:UMI format
#############################################################
### inputs
inDir=$1
i=$2
i_newname=$3
outDir=$4

### set-up
module load samtools/1.9

###
f=${inDir}/${i}

###
echo "starting to work on "${f}
date

samtools view -H ${f} > ${outDir}/tmp.${i}.samheader
samtools view ${f} | sed 's/_/:/' > ${outDir}/tmp.${i}.sambody
cat ${outDir}/tmp.${i}.samheader ${outDir}/tmp.${i}.sambody > ${outDir}/tmp.${i}.sam
echo "done modifying"; date
rm ${outDir}/tmp.${i}.samheader ${outDir}/tmp.${i}.sambody 
samtools view -b $outDir/tmp.${i}.sam | samtools sort -  > $outDir/${i_newname}.bam
echo "done sorting"; date
samtools index -b $outDir/${i_newname}.bam
echo "done indexing"; date
rm $outDir/tmp.${i}.sam
echo "done with "${i}
date
