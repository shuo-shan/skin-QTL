#!/bin/bash

######bsub -Is -q interactive -W 8:00 -n8 -R "span[hosts=1]" /bin/bash
#BSUB -n 1
#BSUB -R rusage[mem=50000]
#BSUB -W 124:00
#BSUB -q long # which queue we want to run in
#BSUB -R span[hosts=1]


module load java/1.8.0_171
module load samtools/1.9
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/dnbseq
###############################################################
id=
cd $Dir/data/bam/barcode_modified
esatPATH=/project/umw_biocore/bin/singleCell/singleCellScripts/esat.v0.1_09.09.16_24.18.umihack.jar
gene=/project/umw_biocore/bin/singleCell/singleCellFiles/hg38_gencode_v34_comprehensive_trans2gene.txt
for i in *$id*.bam; do
  prefix=$(echo $i | cut -d"." -f1)
  java -Xmx40g -jar $esatPATH -alignments $Dir/ESAT/input_names/align_$i -out $prefix -geneMapping $gene -task score3p -wLen 100 -wOlap 50 -wExt 1000 -sigTest .01 -multimap ignore -scPrep
  mv scripture2.log $prefix.ESATscripture2.log
done
###############################################################
#### modify barcode_UMI in read names to barcode:UMI
#id=F49M_IFNg
#cd $Dir/data/bam/softlinks
#for i in *$id*; do
#  f=/nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/DolphinNext_052821/report4643/star/$i
#  outDir=$Dir/data/bam/barcode_modified
#  samtools view -H $f > $outDir/tmp.$i.samheader
#  samtools view $f | sed 's/_/:/' > $outDir/tmp.$i.sambody
#  cat $outDir/tmp.$i.samheader $outDir/tmp.$i.sambody > $outDir/tmp.$i.sam
#  samtools view -b $outDir/tmp.$i.sam | samtools sort - -@ 8 > $outDir/$i
#  samtools index -b -@ 8 $outDir/$i
#  rm $outDir/tmp.$i*
#  echo "done with "$i
#done
#
###############################################################
### split fastq by CS2 barcode
#lib=F27M_F34M
#barcodeMap=$Dir/data/data_info/barcodeMap/$lib
#cd $Dir/data/data_info
#cat samples | grep $lib > this.sample
#cd $Dir/data/fastq/$lib
#while read s;do
#    cd $Dir/data/fastq/$lib/$s
#    for i in $(ls *fq);do
#      gzip $i
#      echo $i
#    done
#    echo ${s}
#done < $Dir/data/data_info/this.sample
#
#lib=F30M_F31M
#barcodeMap=$Dir/data/data_info/barcodeMap/$lib
#cd $Dir/data/data_info
#cat samples | grep $lib > this.sample
#cd $Dir/data/fastq/$lib
#while read s;do
#    cd $Dir/data/fastq/$lib/$s
#    for i in $(ls *fq);do
#      gzip $i
#      echo $i
#    done
#    echo ${s}
#done < $Dir/data/data_info/this.sample
#
#lib=F46M_F47M
#barcodeMap=$Dir/data/data_info/barcodeMap/$lib
#cd $Dir/data/data_info
#cat samples | grep $lib > this.sample
#cd $Dir/data/fastq/$lib
#while read s;do
#    cd $Dir/data/fastq/$lib/$s
#    for i in $(ls *fq);do
#      gzip $i
#      echo $i
#    done
#    echo ${s}
#done < $Dir/data/data_info/this.sample
#
#lib=F49M_F50M
#barcodeMap=$Dir/data/data_info/barcodeMap/$lib
#cd $Dir/data/data_info
#cat samples | grep $lib > this.sample
#cd $Dir/data/fastq/$lib
#while read s;do
#    cd $Dir/data/fastq/$lib/$s
#    for i in $(ls *fq);do
#      gzip $i
#      echo $i
#    done
#    echo ${s}
#done < $Dir/data/data_info/this.sample
#
#lib=F51M_F52M
#barcodeMap=$Dir/data/data_info/barcodeMap/$lib
#cd $Dir/data/data_info
#cat samples | grep $lib > this.sample
#cd $Dir/data/fastq/$lib
#while read s;do
#    cd $Dir/data/fastq/$lib/$s
#    for i in $(ls *fq);do
#      gzip $i
#      echo $i
#    done
#    echo ${s}
#done < $Dir/data/data_info/this.sample
#
#lib=F55M_F56M
#barcodeMap=$Dir/data/data_info/barcodeMap/$lib
#cd $Dir/data/data_info
#cat samples | grep $lib > this.sample
#cd $Dir/data/fastq/$lib
#while read s;do
#    cd $Dir/data/fastq/$lib/$s
#    for i in $(ls *fq);do
#      gzip $i
#      echo $i
#    done
#    echo ${s}
#done < $Dir/data/data_info/this.sample
