### 05/25/2021
### script for CELseq2 data processing
### working in /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/scripts
### set-up
bsub -Is -q interactive -W 8:00 -n1 -R rusage[mem=450000] -R "span[hosts=1]" /bin/bash
module load java/1.8.0_171
module load condas/2018-05-11
source activate sshan_isoform
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/dnbseq

##############################################################
### split fastq by CS2 barcode
lib=F27M_F34M
barcodeMap=$Dir/data/data_info/barcodeMap/$lib
cd $Dir/data/data_info
cat samples | grep $lib > this.sample
cd $Dir/data/fastq/$lib
while read s;do
    gunzip ${s}_1.fq.gz && gunzip ${s}_2.fq.gz
    java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_1.fq F2=${s}_2.fq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
    gzip ${s}_1.fq && gzip ${s}_2.fq
    cd $Dir/data/fastq/$lib/$s
    for i in $(ls */*fq);do
      gzip $i
      echo $i
    done
    echo ${s}
done < $Dir/data/data_info/this.sample
##############################################################
### run on DolphinNext pipeline
# create soft links of all fq.gz files in /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/data/fastq/softlink
# map to hg38
# paired-end
# run STAR, FeatureCounts_after_STAR, RSEM, IGV_IDF_conversion, RSeQC, Quality_Filtering, FastQC
# https://dolphinnext.umassmed.edu/index.php?np=3&id=4628
# follow-up: didn't work because P1 reads are only 3 bp long!
# run pipeline again for just P2: /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/DolphinNext_052821
# https://dolphinnext.umassmed.edu/index.php?np=3&id=4643

##############################################################
### run star-aligned bam through ESAT
# /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/DolphinNext_052821/report4643/star/*.bam
# create soft links of all bam files in /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/data/bam/softlink
# change read name to have :barcode:UMI format
cd $Dir/data/bam/softlinks
for i in *; do
  f=/nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/DolphinNext_052821/report4643/star/$i
  outDir=$Dir/data/bam/barcode_modified
  samtools view -H $f > $outDir/tmp.$i.samheader
  samtools view $f | sed 's/_/:/' > $outDir/tmp.$i.sambody
  cat $outDir/tmp.$i.samheader $outDir/tmp.$i.sambody > $outDir/tmp.$i.sam
  samtools view -b $outDir/tmp.$i.sam | samtools sort - -@ 8 > $outDir/$i
  samtools index -b -@ 8 $outDir/$i
  rm $outDir/tmp.$i*
  echo "done with "$i
done
# create file name and path for ESAT input
cd $Dir/data/bam/barcode_modified
for i in *_sorted.bam; do
  id=$(echo $i | sed "s/_sorted.bam//g" | cut -d"_" -f2-3)
  echo -e "$id"'\t'/nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/data/bam/barcode_modified/"$i" > $Dir/ESAT/input_names/align_$i
done
# run ESAT
cd $Dir/data/bam/barcode_modified
esatPATH=/project/umw_biocore/bin/singleCell/singleCellScripts/esat.v0.1_09.09.16_24.18.umihack.jar
gene=/project/umw_biocore/bin/singleCell/singleCellFiles/hg38_gencode_v34_comprehensive_trans2gene.txt
for i in *.bam; do
  prefix=$(echo $i | cut -d"." -f1)
  java -Xmx40g -jar $esatPATH -alignments $Dir/ESAT/input_names/align_$i -out $prefix -geneMapping $gene -task score3p -wLen 100 -wOlap 50 -wExt 1000 -sigTest .01 -multimap ignore -scPrep 
  mv scripture2.log $prefix.ESATscripture2.log
done
mv *.gene.txt /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/ESAT/output/gene
mv *.window.txt /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/ESAT/output/window
mv *ESATscripture2.log /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/ESAT/scripture
mv *umi.distributions.txt /nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/ESAT/umi_distributions
# organize ESAT output
cd $Dir/ESAT/output/gene
# paste is appropriate here b/c all files have the same row names
paste * > temp.MEL.gene.txt
# get rid of duplicated gene symbol, chr, strand columns from pasting:
keep=$(head -1 temp.MEL.gene.txt | tr '\t' '\n' | grep -v -n -E 'Symbol|chr|strand' | cut -d":" -f1 | tr '\n' ',' | sed 's/\(.*\),/\1/' | sed 's/4,/1,2,3,4,/')
# get rid of sample barcode from sample name
cat temp.MEL.gene.txt | cut -f$keep | sed '1 s/:[ATCG]\{6\}//g' > MEL.gene.txt
rm temp.MEL.gene.txt
chmod 777 MEL.gene.txt


















