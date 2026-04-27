#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=120000]
#BSUB -q long
#BSUB -W 21:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2/%J%I.AmazonBackup.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2/%J%I.AmazonBackup.err"
module load bcl2fastq/2.17.1.14
# 08_12_21 sequencing data with Ken, Wei
# Illumina NextSeq550, pair-end, single-index 8bp for ChIPseq and 6bp for Celseq2 RPI barcode
#indir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2/basespace
#samplesheet1=$indir/SampleSheet_ChIPseq.csv
#samplesheet2=$indir/SampleSheet_Celseq2.csv
#outdir=$indir/fastqs
#bcl2fastq \
#--runfolder-dir $indir \
#--sample-sheet $samplesheet1 \
#--output-dir $outdir \
#--interop-dir $outdir \
#--use-bases-mask Y*,I8n*,Y* \
#--mask-short-adapter-reads 0 \
#--minimum-trimmed-read-length 0 \
#--barcode-mismatches 1
#echo "done with bcl2fastq for ChIPseq data"
#date
#bcl2fastq \
#--runfolder-dir $indir \
#--sample-sheet $samplesheet2 \
#--output-dir $outdir \
#--interop-dir $outdir \
#--use-bases-mask Y*,I6n*,Y* \
#--mask-short-adapter-reads 0 \
#--minimum-trimmed-read-length 0 \
#--barcode-mismatches 1
#echo "done with bcl2fastq for Celseq2 data"
#date
# resulted in:
#17128170.bcl2fq.err and 17128170.bcl2fq.out
##########################################
### Note to self:
### Prior to this bcl2fastq step, all data were downloaded to cluster from galaxy
#     Step 1: in galaxy, run basemount --unmount /home/shans/BaseSpace, then basemount /home/shans/BaseSpace
#     Step 2: go to /home/shans/BaseSpace/Runs/, then go to my run folder, go to Files
#     Step 3: rsync -av *.xml *.txt to_path_in_cluster/
#             rsync -av Data/ to_path_in_cluster/Data/
### Back in cluster
#     Step 1: go to to_path_in_cluster/fastqs folder, modify samplesheet
#     Step 2: modify this bcl2fastq file
##########################################
### merge
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2/basespace/fastqs/chip_rnaseq_08-12-21
#cd $dir
##ls | cut -d"_" -f1-3 | sort | uniq | tr '\n' ' '
#for x in Celseq2_1_S1 Celseq2_2_S2 Celseq2_3_S3 Celseq2_4_S4 Celseq2_5_S5 Celseq2_6_S6 ChIP_1_S1 ChIP_2_S2 ChIP_3_S3 ChIP_4_S4 ChIP_5_S5 ChIP_6_S6 ChIP_7_S7;do
#  #zcat ${dir}/${x}_L001_R1_001.fastq.gz ${dir}/${x}_L002_R1_001.fastq.gz ${dir}/${x}_L003_R1_001.fastq.gz ${dir}/${x}_L004_R1_001.fastq.gz > $dir/../laneMerged/${x}_R1.fastq
#  #gzip $dir/../laneMerged/${x}_R1.fastq
#  #rm $dir/../laneMerged/${x}_R1.fastq
#  zcat ${dir}/${x}_L001_R2_001.fastq.gz ${dir}/${x}_L002_R2_001.fastq.gz ${dir}/${x}_L003_R2_001.fastq.gz ${dir}/${x}_L004_R2_001.fastq.gz > $dir/../laneMerged/${x}_R2.fastq
#  gzip $dir/../laneMerged/${x}_R2.fastq
#  rm $dir/../laneMerged/${x}_R2.fastq
#  echo "done with "$lib
#done
#########################################
## split library into samples in Celseq2 by CS2 barcode
module load java/1.8.0_171
module load condas/2018-05-11
source activate sshan_isoform
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2
for s in Celseq2_1_S1 Celseq2_2_S2 Celseq2_3_S3 Celseq2_4_S4 Celseq2_5_S5 Celseq2_6_S6;do
  echo "working on" ${s}
  date
  cd $Dir/basespace/fastqs/laneMerged
  barcodeMap=$Dir/data_info_Celseq2/barcodeMap/${s}
  gunzip ${s}_R1.fastq.gz && gunzip ${s}_R2.fastq.gz
  java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_R1.fastq F2=${s}_R2.fastq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
  gzip ${s}_R1.fastq && gzip ${s}_R2.fastq
  cd $Dir/basespace/fastqs/laneMerged/${s}
  for i in $(ls *fq);do
    gzip ${i}
    echo ${i}
    date
  done
  echo ${s}
  date
done

### check md5sum then backup on Amazon
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2
cd $Dir/basespace/fastqs/laneMerged
for i in `find . -type f`; do md5sum $i > $i.md5sum; done
cd ..
/project/umw_biocore/bin/amazonBackup.bash laneMerged s3://biocorebackup/garberlab/human/skin/eQTL/Nextseq550/20210812_ChIPseq_Celseq2/
