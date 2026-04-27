#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=120000]
#BSUB -q long
#BSUB -W 21:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210903_ATACseq_Celseq2/%J%I.bcl2fq_laneMerge.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210903_ATACseq/%J%I.bcl2fq_laneMerge.err"
module load bcl2fastq/2.17.1.14
# 09_03_21 sequencing data with Ken
# Illumina NextSeq550, pair-end, single-index 8bp for ATACseq, Read1(80bp), Read2(80bp)
indir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210903_ATACseq/basespace
samplesheet=$indir/SampleSheet.csv
outdir=$indir/fastqs
bcl2fastq \
--runfolder-dir $indir \
--sample-sheet $samplesheet \
--output-dir $outdir \
--interop-dir $outdir \
--use-bases-mask Y*,I8n*,Y* \
--mask-short-adapter-reads 0 \
--minimum-trimmed-read-length 0 \
--barcode-mismatches 1
echo "done with bcl2fastq for ATACseq data"
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
#########################################
## merge
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210903_ATACseq/basespace/fastqs/atac_rnaseq_09-3-21
#cd $dir
##ls | cut -d"_" -f1-3 | sort | uniq | tr '\n' ' '
#for x in ATACseq_1_S1 ATACseq_2_S2;do
#  zcat ${dir}/${x}_L001_R1_001.fastq.gz ${dir}/${x}_L002_R1_001.fastq.gz ${dir}/${x}_L003_R1_001.fastq.gz ${dir}/${x}_L004_R1_001.fastq.gz > $dir/../laneMerged/${x}_R1.fastq
#  gzip $dir/../laneMerged/${x}_R1.fastq
#  rm $dir/../laneMerged/${x}_R1.fastq
#  zcat ${dir}/${x}_L001_R2_001.fastq.gz ${dir}/${x}_L002_R2_001.fastq.gz ${dir}/${x}_L003_R2_001.fastq.gz ${dir}/${x}_L004_R2_001.fastq.gz > $dir/../laneMerged/${x}_R2.fastq
#  gzip $dir/../laneMerged/${x}_R2.fastq
#  rm $dir/../laneMerged/${x}_R2.fastq
#  echo "done with "$lib
#done

#### check md5sum then backup on Amazon
#Dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2
#cd $Dir/basespace/fastqs/laneMerged
#for i in `find . -type f`; do md5sum $i > $i.md5sum; done
#cd ..
#/project/umw_biocore/bin/amazonBackup.bash laneMerged s3://biocorebackup/garberlab/human/skin/eQTL/Nextseq550/20210812_ChIPseq_Celseq2/
