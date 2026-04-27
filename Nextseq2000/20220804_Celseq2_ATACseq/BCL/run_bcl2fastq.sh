#!/bin/bash
#BSUB -n 5
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 24:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220804_Celseq2_ATACseq/BCL/bcl2fq_ATACseq_%J%I.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220804_Celseq2_ATACseq/BCL/bcl2fq_ATACseq_%J%I.err"

module load bcl2fastq2/2.20.0

# 08_04_22 Nextseq2000 sequencing data prepared by Crystal Shan
# Illumina NextSeq2k, pair-end, single-index 8bp
# first job for Celseq2 RNAseq libraries (6bp index)
#indir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220804_Celseq2_ATACseq/BCL
#samplesheet=$indir/SampleSheet_Celseq2.csv
#outdir=$indir/../fastqs_RNAseq
#
#bcl2fastq \
#--runfolder-dir $indir \
#--output-dir $outdir \
#--interop-dir $outdir \
#--stats-dir $outdir \
#--sample-sheet $samplesheet \
#--use-bases-mask Y*,I6n*,Y* \
#--mask-short-adapter-reads 0 \
#--minimum-trimmed-read-length 0 \
#--barcode-mismatches 1 \
#--loading-threads 1 \
#--processing-threads 4 \
#--writing-threads 1

# another run for ATACseq libraries (8bp index)
indir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220804_Celseq2_ATACseq/BCL
samplesheet=$indir/SampleSheet_ATAC.csv
outdir=$indir/../fastqs_ATACseq

bcl2fastq \
--runfolder-dir $indir \
--output-dir $outdir \
--interop-dir $outdir \
--stats-dir $outdir \
--sample-sheet $samplesheet \
--use-bases-mask Y*,I8n*,Y* \
--mask-short-adapter-reads 0 \
--minimum-trimmed-read-length 0 \
--barcode-mismatches 1 \
--loading-threads 1 \
--processing-threads 4 \
--writing-threads 1


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

### Nextseq200 on basespace uses BCL convert, which doesn't provide fastq files of undetermined indices
#     Back in cluster, if I want to manually paarse BCL to fastq, uses bcl2fastq2, bcl2fastq2/2.20.0, and the samplesheet
#     SHALL NOT be downloaded from basespace (because it's BCLconvert format), I should use the same samplesheet
#     format as Nextseq550.
