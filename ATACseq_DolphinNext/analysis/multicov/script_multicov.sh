#!/bin/bash
#BSUB -n 26
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=18000]
#BSUB -q long
#BSUB -W 121:00
### comment

#### packages & env
module load bedtools/2.29.2

### 1. merge all narrowPeaks
wd=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis
pth=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/narrowPeak
cat ${pth}/*.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n > ${wd}/temp.one
bedtools merge -i ${wd}/temp.one > ${wd}/merged.bed
rm ${wd}/temp.one

### 2. perform multicov on merged bed
# use function_calc_cov.sh
while read c; do
	echo ${c} | bsub -W 8:00 -n 9 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long -o "/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/multicov/%J%I.out" -e "/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/multicov/%J%I.err"
done < commands.txt

### 3. paste all bedgraphs together
wd=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/multicov
cd ${wd}
# get the first 3 columns of peaks
bedF=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/merged.bed
printf "chr\tstart\tend\n" >> temp.header
cat temp.header ${bedF} > temp.bed
# paste to master file
ls ${wd}/bedgraph | sed 's/.merged.bg//g' > temp.names
while read line;do
  f=${wd}/bedgraph/${line}.merged.bg
  echo ${line} > temp.${line}.this
  cat ${f} | cut -f4 >> temp.${line}.this
  echo "done with "$line
done < temp.names
paste temp.bed temp*this > multicov.merged.txt
rm temp*

