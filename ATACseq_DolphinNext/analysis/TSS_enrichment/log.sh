#!/bin/bash
# Transcription Start Site (TSS) Enrichment Score. 
# shuo.shan@umassmed.edu
# Nov2021

#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=4080]
#BSUB -q long
#BSUB -W 8:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/TSS_enrichment/%J%I.bin.out" 
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/TSS_enrichment/%J%I.bin.err" 
### comment

#### user input
sample=$1
#sample=ATAC_F25M_PBS_S1
sampleF=${bamDir}/${sample}.bam
#### packages & env
module load bedtools/2.29.2
module load samtools/1.9
wd=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/TSS_enrichment
bamDir=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/bam
#### 

# 1. get start 1bp of every gene (TSS)
#geneF=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/gencode_v34/genes/genes.bed
#genomeF=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
#cat ${geneF} | awk '{OFS="\t"}{print $1,$2,$2+1}' > TSS.bed
#cat ${geneF} | awk '{OFS="\t"}{if ($2>2000) print $1,$2-2000,$2+2000}' > TSSflank
#cat TSSflank | sed 's/\t/_/g' > TSSflank.name
#paste TSSflank TSSflank.name > TSSflank.bed
#rm TSSflank.name TSSflank 

# 2. get 100bp windows per TSS flanking region.
#while read entry;do
#  name=$(echo ${entry} | cut -d" " -f4)
#  yes ${name} | head -n 40 > temp.${name}.txt
#  echo ${entry} | tr " " "\t" > temp.${name}.bed
#  bedtools makewindows -b temp.${name}.bed -w 100 > temp.${name}.binned
#  paste temp.${name}.binned temp.${name}.txt | awk '{OFS="\t"}{print $1,$2,$3,$4"_"NR}' >> TSS.binned.bed
#  echo ${name} >> finished.txt
#  rm temp.${name}.txt temp.${name}.bed temp.${name}.binned
#done < TSSflank.bed

# 3. get coverage for every window
#cd ${bamDir}
#samtools index ${sampleF}
#cd ${wd}
#bedtools multicov -bams ${sampleF} -bed TSS.binned.bed > TSS.binned.cov.${sample}.bg
while read c; do
  tag=$(echo ${c} | cut -d" " -f3)
  echo ${c} | bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long -o "./%J%I.coverage.${tag}.out"  -e "./%J%I.coverage.${tag}.err" -J ${tag}
done < commands.txt
