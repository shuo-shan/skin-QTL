#!/bin/bash
# merge peaks across all samples
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16000]"
#BSUB -W 08:00
#BSUB -q long
#BSUB -J perSampleCov[1-36]
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/log/countCoveragePeaks_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/log/countCoveragePeaks_QTL_%J_%I.err"

module load bedtools/2.30.0
module load samtools/1.16.1

# --------------------- #
# Perform per-sample coverage counting
# --------------------- #
MASTER=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/master_peaks.bed
BAM_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped
COL_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/cols
mkdir -p ${COL_DIR}

# Build sorted BAM array (same order every time)
BAM_LIST=($(ls ${BAM_DIR}/*.nodup.bam | sort))

# LSF array index is 1-based
IDX=$((LSB_JOBINDEX - 1))
BAM="${BAM_LIST[$IDX]}"
SAMPLE=$(basename ${BAM} .nodup.bam)

echo "Processing sample ${LSB_JOBINDEX}: ${SAMPLE}"

# Extract only the count column (col 4 of bedtools coverage output)
bedtools coverage \
    -sorted \
    -counts \
    -g /share/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.chrom.sizes \
    -a ${MASTER} \
    -b ${BAM} \
    | awk '{print $NF}' \
    > ${COL_DIR}/${SAMPLE}.counts.txt

echo "Done: ${SAMPLE}"; date
