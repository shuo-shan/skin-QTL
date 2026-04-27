#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=81920]
#BSUB -q long
#BSUB -W 21:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210624_ATAC_Celseq2/%J%I.zip.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210624_ATAC_Celseq2/%J%I.zip.err"

module load bcl2fastq/2.17.1.14

# 06_24_21 sequencing data with Ken, Wei
# Illumina NextSeq550, pair-end, single-index 8bp for ATACseq and 6bp for Celseq2 RPI barcode
indir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210624_ATAC_Celseq2
samplesheet1=$indir/SampleSheet_ATAC.csv
samplesheet2=$indir/SampleSheet_Celseq2.csv
outdir=$indir/fastqs

#bcl2fastq \
#--runfolder-dir $indir \
#--sample-sheet $samplesheet1 \
#--output-dir $outdir \
#--interop-dir $outdir \
#--use-bases-mask Y*,I8n*,Y* \
#--mask-short-adapter-reads 0 \
#--minimum-trimmed-read-length 0 \
#--barcode-mismatches 1
#
#
#bcl2fastq \
#--runfolder-dir $indir \
#--sample-sheet $samplesheet2 \
#--output-dir $outdir \
#--interop-dir $outdir \
#--use-bases-mask Y*,I6n*,Y* \
#--mask-short-adapter-reads 0 \
#--minimum-trimmed-read-length 0 \
#--barcode-mismatches 1
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
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210624_ATAC_Celseq2/fastqs
#
#cd $dir/skin_eQTL_ATAC
#for x in ATAC_F25F_IFN_S2 ATAC_F25F_PBS_S1 ATAC_F49F_IFN_S4 ATAC_F49F_PBS_S3 ATAC_F55F_IFN_S6 ATAC_F55F_PBS_S5;do
#  lib=$(echo ${x} | cut -d"_" -f1,2,3)
#  sample=$(echo $x | cut -d"_" -f2,3,4)
#  echo "lib is " $lib
#  zcat $dir/skin_eQTL_ATAC/${lib}/${sample}_L001_R1_001.fastq.gz $dir/$lib/${sample}_L002_R1_001.fastq.gz $dir/$lib/${sample}_L003_R1_001.fastq.gz $dir/$lib/${sample}_L004_R1_001.fastq.gz > $dir/merged/${lib}_R1.fastq.gz
#  zcat $dir/skin_eQTL_ATAC/${lib}/${sample}_L001_R2_001.fastq.gz $dir/$lib/${sample}_L002_R2_001.fastq.gz $dir/$lib/${sample}_L003_R2_001.fastq.gz $dir/$lib/${sample}_L004_R2_001.fastq.gz > $dir/merged/${lib}_R2.fastq.gz
#  echo "done with "$lib
#done
#
#cd $dir/skin_eQTL_CelSeq2
#for x in Celseq2_F61K_F62K_S1 Celseq2_F62M_F63M_S3 Celseq2_F63K_F22F_S2;do
#  lib=$(echo ${x} | cut -d"_" -f1,2,3)
#  sample=$(echo $x | cut -d"_" -f2,3,4)
#  echo "lib is " $lib
#  zcat $dir/skin_eQTL_CelSeq2/${lib}/${sample}_L001_R1_001.fastq.gz $dir/$lib/${sample}_L002_R1_001.fastq.gz $dir/$lib/${sample}_L003_R1_001.fastq.gz $dir/$lib/${sample}_L004_R1_001.fastq.gz > $dir/merged/${lib}_R1.fastq.gz
#  zcat $dir/skin_eQTL_CelSeq2/${lib}/${sample}_L001_R2_001.fastq.gz $dir/$lib/${sample}_L002_R2_001.fastq.gz $dir/$lib/${sample}_L003_R2_001.fastq.gz $dir/$lib/${sample}_L004_R2_001.fastq.gz > $dir/merged/${lib}_R2.fastq.gz
#  echo "done with "$lib
#done
#

#### check md5sum then backup on Amazon
#cd $dir/merged
#for i in `find . -type f`; do md5sum $i > $i.md5sum; done
#cd ..
#/project/umw_biocore/bin/amazonBackup.bash merged s3://biocorebackup/garberlab/human/skin/eQTL/Nextseq550/20210624_ATAC_Celseq2

##########################################
### split library into samples in Celseq2 by CS2 barcode
# Celseq2_F61K_F62K Celseq2_F62M_F63M Celseq2_F63K_F22F
module load java/1.8.0_171
module load condas/2018-05-11
source activate sshan_isoform
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210624_ATAC_Celseq2
for s in Celseq2_F61K_F62K Celseq2_F62M_F63M Celseq2_F63K_F22F;do
#  cd $Dir/fastqs/merged
#  barcodeMap=$Dir/data_info_Celseq2/barcodeMap/${s}
#  gunzip ${s}_R1.fastq.gz && gunzip ${s}_R2.fastq.gz
#  java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_R1.fastq F2=${s}_R2.fastq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
#  gzip ${s}_R1.fastq && gzip ${s}_R2.fastq
  cd $Dir/fastqs/merged/${s}
  for i in $(ls *fq);do
    gzip ${i}
    echo ${i}
  done
  echo ${s}
done


