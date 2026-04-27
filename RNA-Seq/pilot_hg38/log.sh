#!/bin/bash
### 07/14/2021

#BSUB -n 1 
#BSUB -R rusage[mem=120000]
#BSUB -W 124:00
#BSUB -q long
#BSUB -R span[hosts=1]
##bsub -Is -q interactive -W 8:00 -n1 -R rusage[mem=450000] -R "span[hosts=1]" /bin/bash

### script for CELseq2 data processing
### working in /nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/pilot_hg38
### set-up
module load java/1.8.0_171
module load samtools/1.9
module load condas/2018-05-11
source activate sshan_isoform
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/pilot_hg38
##############################################################
# These files are from the pilot study done on F25 samples for all 3 celltypes.
# The goal is to identify the appropriate IFNg stimulation duration and amount.
# It was orignally mapped to hg19 genome, so I am re-mapping the fastq to hg38.
# 
##############################################################
#### merge fastq files
#for bc in AGACTC CATGAG CAGATC TCACAG GTCTAG GTTGCA ACCATG;do
#  for ct in MEL FRB KRT;do
#    for mate in p1 p2;do
#      zcat /nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/fastqs/RNASeq/celseq_${ct}*/${bc}*${mate}.fq.gz | gzip > ${Dir}/fastq/celseq_${ct}_${bc}.${mate}.fq.gz
#      echo done with celseq_${ct}_${bc}.${mate}.fq.gz
#    done 
#  done
#done
#
##############################################################
### run on DolphinNext pipeline
# create soft links of all fq.gz files in /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/data/fastq/softlink
# map to hg38 (ver28)
# single-end with only p2 reads. RNAseq pipeline ver5.
# run STAR, FeatureCounts_after_STAR, RSEM, IGV_IDF_conversion, RSeQC, Quality_Filtering, FastQC
# https://dolphinnext.umassmed.edu/index.php?np=3&id=4907


##############################################################
### run star-aligned bam through ESAT
# /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/DolphinNext_052821/report4643/star/*.bam
# create soft links of all bam files in /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/data/bam/softlink
# change read name to have :barcode:UMI format
#cd ${Dir}/DolphinNext/report4907/star
#for i in *.bam; do
#  f=${Dir}/DolphinNext/report4907/star/${i}
#  outDir=$Dir/bam_bcmodified
#  samtools view -H $f > $outDir/tmp.$i.samheader
#  samtools view $f | sed 's/_/:/' > $outDir/tmp.$i.sambody
#  cat $outDir/tmp.$i.samheader $outDir/tmp.$i.sambody > $outDir/tmp.$i.sam
#  samtools view -b $outDir/tmp.$i.sam | samtools sort - -@ 8 > $outDir/$i
#  samtools index -b -@ 8 $outDir/$i
#  rm $outDir/tmp.$i*
#  echo "done with "$i
#done
## create file name and path for ESAT input
#cd ${Dir}/bam_bcmodified
#for i in *_sorted.bam; do
#  id=$(echo $i | sed "s/_sorted.bam//g" | cut -d"_" -f2-3)
#  echo -e "F25_${id}"'\t'/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/pilot_hg38/bam_bcmodified/"$i" > ${Dir}/ESAT/input_names/align_${i}
#done
# run ESAT
cd ${Dir}/bam_bcmodified
esatPATH=/project/umw_biocore/bin/singleCell/singleCellScripts/esat.v0.1_09.09.16_24.18.umihack.jar
gene=/project/umw_biocore/bin/singleCell/singleCellFiles/hg38_gencode_v34_comprehensive_trans2gene.txt
for i in *.bam; do
  prefix=$(echo $i | cut -d"." -f1)
  java -Xmx100g -jar $esatPATH -alignments $Dir/ESAT/input_names/align_$i -out $prefix -geneMapping $gene -task score3p -wLen 100 -wOlap 50 -wExt 1000 -sigTest .01 -multimap ignore -scPrep 
  mv scripture2.log $prefix.ESATscripture2.log
  echo "done with "$i
done
mv *.gene.txt ${Dir}/ESAT/output/gene
mv *.window.txt ${Dir}/ESAT/output/window
mv *ESATscripture2.log ${Dir}/ESAT/scripture
mv *umi.distributions.txt ${Dir}/ESAT/umi_distributions
## organize ESAT output
cd $Dir/ESAT/output/gene
# paste is appropriate here b/c all files have the same row names
bcMap=/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/fastqs/RNASeq/barcodeMap
cd ${Dir}/ESAT/output/gene
for bc in AGACTC CATGAG CAGATC TCACAG GTCTAG GTTGCA ACCATG;do
  for ct in MEL FRB KRT;do
    condition=$(grep ${bc} ${bcMap} | cut -f2)
    cp celseq_${ct}_${bc}_sorted.gene.txt celseq_${ct}_${condition}_sorted.gene.txt
    echo "done with "${ct} ${bc}
  done
done
## get rid of duplicated gene symbol, chr, strand columns from pasting:
paste *${ct}*${bc}* > temp.${ct}.${bc}.gene.txt
keep=$(head -1 temp.MEL.gene.txt | tr '\t' '\n' | grep -v -n -E 'Symbol|chr|strand' | cut -d":" -f1 | tr '\n' ',' | sed 's/\(.*\),/\1/' | sed 's/4,/1,2,3,4,/')
## get rid of sample barcode from sample name
cat temp.MEL.gene.txt | cut -f$keep | sed '1 s/:[ATCG]\{6\}//g' > MEL.gene.txt
#rm temp.MEL.gene.txt
#chmod 777 MEL.gene.txt
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#


