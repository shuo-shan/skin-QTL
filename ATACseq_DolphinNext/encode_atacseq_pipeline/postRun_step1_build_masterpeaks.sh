#!/bin/bash
# merge peaks across all samples
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=100000]"
#BSUB -W 08:00
#BSUB -q long
#BSUB -J makeMasterPeak
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/log/mergePeaks_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/log/mergePeaks_QTL_%J_%I.err"

module load samtools/1.16.1
module load macs2/2.2.7.1
module load bedtools/2.30.0

# --------------------- #
# Make summit peaks
# --------------------- #

MACS2_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2
OUT_DIR=${MACS2_DIR}/summit_peaks
mkdir -p ${OUT_DIR}

EXTEND=250  # ±250bp around summit = 500bp peak
CHROM_SIZES=/share/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.chrom.sizes

for summit_file in ${MACS2_DIR}/*_summits.bed; do
    sample=$(basename ${summit_file} _summits.bed)
    echo "Processing ${sample}..."

    # Extend summit ± 250bp, clip to chrom boundaries
    bedtools slop -i ${summit_file} \
        -g ${CHROM_SIZES} \
        -b ${EXTEND} \
        > ${OUT_DIR}/${sample}.summit250.bed

    echo "Done: ${sample}"
done

echo "All summit peaks generated"; date

# --------------------- #
# Build the master peak list
# --------------------- #
SUMMIT_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2/summit_peaks
OUT_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks
mkdir -p ${OUT_DIR}

# Concatenate all summit-centered peaks, sort, merge overlapping ones
cat ${SUMMIT_DIR}/*.summit250.bed | grep -v "^chrM" | grep -v "^chrUn" | grep -v "_random" | sort -k1,1V -k2,2n |\
	bedtools merge -i - > ${OUT_DIR}/master_peaks_raw.bed

# Add peak names (ATACseq_peak_1, ATACseq_peak_2, ...)
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "ATACseq_peak_"NR}' \
    ${OUT_DIR}/master_peaks_raw.bed \
    > ${OUT_DIR}/master_peaks_prefiltered.bed

echo "Total peaks in master list:"
wc -l ${OUT_DIR}/master_peaks_prefiltered.bed

# Remove blacklisted regions
BLACKLIST=/pi/manuel.garber-umw/human/skin/eQTLs/literature/blacklist/hg38.blacklist.bed
bedtools intersect -v \
    -a ${OUT_DIR}/master_peaks_prefiltered.bed \
    -b ${BLACKLIST} \
    > ${OUT_DIR}/master_peaks_filtered.bed

# Rename so peak IDs are clean and sequential AFTER filtering
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "ATACseq_peak_"NR}' \
    ${OUT_DIR}/master_peaks_filtered.bed \
    > ${OUT_DIR}/master_peaks.bed

echo "Total peaks in master list after blacklist filtering:"
wc -l ${OUT_DIR}/master_peaks.bed

rm ${OUT_DIR}/master_peaks_raw.bed ${OUT_DIR}/master_peaks_prefiltered.bed ${OUT_DIR}/master_peaks_filtered.bed



echo "Done building master peak list"; date


