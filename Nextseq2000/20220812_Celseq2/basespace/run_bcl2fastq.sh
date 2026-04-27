#!/bin/bash
#BSUB -n 9
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 8:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220812_Celseq2/basespace/bcl2fq_%J%I.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220812_Celseq2/basespace/bcl2fq_%J%I.err"

module load bcl2fastq2/2.20.0

# 08_12_22 Nextseq2000 sequencing data prepared by Crystal Shan
# Illumina NextSeq2k, pair-end, single-index 6bp
indir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220812_Celseq2/basespace
samplesheet=$indir/SampleSheet_Celseq2.csv
outdir=$indir/../fastq_RNAseq

bcl2fastq \
--runfolder-dir $indir \
--output-dir $outdir \
--interop-dir $outdir \
--stats-dir $outdir \
--sample-sheet $samplesheet \
--use-bases-mask Y*,I6n*,Y* \
--mask-short-adapter-reads 0 \
--minimum-trimmed-read-length 0 \
--barcode-mismatches 1 \
--loading-threads 8 \
--processing-threads 8 \
--writing-threads 8


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
