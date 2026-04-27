#!/bin/bash
# merge peaks across all samples
# call Rscript to normalize count, define open peaks, define differential accessibility peaks, compile big table
#BSUB -n 3
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=60000]"
#BSUB -W 08:00
#BSUB -q long
#BSUB -J ATACseqDE
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/log/ATACseqDAC_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/log/ATACseqDAC_%J_%I.err"


module load bedtools/2.30.0

dir=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/
cd ${dir}

MACS2_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2
MASTER=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/master_peaks.bed
OUT_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/peak_overlap
mkdir -p ${OUT_DIR}

# Use the _peaks.narrowPeak files (unzipped, already there)
NARROW_LIST=($(ls ${MACS2_DIR}/*.clipped.narrowPeak.gz | sort))

echo "Found ${#NARROW_LIST[@]} narrowPeak files"; date

# For each sample, get binary overlap with master peaks (1 = overlaps, 0 = doesn't)
for narrow in "${NARROW_LIST[@]}"; do
    sample=$(basename ${narrow} .clipped.narrowPeak.gz)
    echo "Processing ${sample}..."; date

    bedtools intersect \
        -a ${MASTER} \
        -b ${narrow} \
        -c \
        | awk '{print ($NF > 0) ? 1 : 0}' \
        > ${OUT_DIR}/${sample}.overlap.txt

    echo "Done: ${sample}"
done

echo "All overlaps done"; date

# Build header
HEADER="peak_name"
for narrow in "${NARROW_LIST[@]}"; do
    sample=$(basename ${narrow} .clipped.narrowPeak.gz)
    HEADER="${HEADER}\t${sample}"
done

# Paste all columns together with peak names
# Extract peak names from master
cut -f4 ${MASTER} > ${OUT_DIR}/peak_names.txt

# Collect overlap columns in same sorted order
COL_FILES=()
for narrow in "${NARROW_LIST[@]}"; do
    sample=$(basename ${narrow} .clipped.narrowPeak.gz)
    COL_FILES+=("${OUT_DIR}/${sample}.overlap.txt")
done

echo -e "${HEADER}" > ${OUT_DIR}/master_peaks_overlap.txt
paste ${OUT_DIR}/peak_names.txt "${COL_FILES[@]}" >> ${OUT_DIR}/master_peaks_overlap.txt

echo "Done! Output: ${OUT_DIR}/master_peaks_overlap.txt"; date
echo "Dimensions:"
wc -l ${OUT_DIR}/master_peaks_overlap.txt
head -2 ${OUT_DIR}/master_peaks_overlap.txt

# ------------------------------------------------------------------------
# ------- Performing ATACseq now ----------
# ------------------------------------------------------------------------


echo "performing Rscript for ATACseq merged count matrix now"
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${dir}/postRun_step4_compile_table_for_accessibility_and_dynamics.R

