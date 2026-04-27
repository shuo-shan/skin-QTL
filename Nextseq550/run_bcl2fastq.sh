#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=8192]
#BSUB -q long
#BSUB -W 8:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210608_ATAC_ChIPseq/%J%I.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210608_ATAC_ChIPseq/%J%I.err"

module load bcl2fastq/2.17.1.14

# 06_08_21 sequencing data with Ken, Wei
# Illumina NextSeq550, pair-end, single-index 8bp
indir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210608_ATAC_ChIPseq
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

