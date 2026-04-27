#!/bin/bash
# this script in ENCODE ATACseq pipeline performs post-alignment filtering, removes PCR duplicates, and calculates library complexity.

# -------------------------------------------------------------------------
# set-up
prefix=$1
dir=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline
filtered_alignment_dir=${dir}/filtered_alignment
cd ${dir}

source activate ATACseq
module load samtools
module load bedtools
module load openjdk

# =============
# Mark duplicates
# =============
echo "marking duplicates for ${prefix}"; date
FILT_BAM_PREFIX=${prefix}.filt
FILT_BAM_FILE=${filtered_alignment_dir}/${prefix}.filt.bam
TMP_FILT_BAM_FILE=${filtered_alignment_dir}/${prefix}.filt.dupmark.bam
DUP_FILE_QC=${filtered_alignment_dir}/${prefix}.filt.dup.qc
MARKDUP="/pi/manuel.garber-umw/sshan/scripts/picard.jar" # picard 2.27.4

ls -lh ${FILT_BAM_FILE}
java -Xmx4G -jar ${MARKDUP} MarkDuplicates -INPUT ${FILT_BAM_FILE} -OUTPUT ${TMP_FILT_BAM_FILE} -METRICS_FILE ${DUP_FILE_QC} -VALIDATION_STRINGENCY LENIENT -ASSUME_SORTED true -REMOVE_DUPLICATES false -READ_NAME_REGEX "null"

mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}

# ============================
# Remove duplicates
# Index final position sorted BAM
# Create final name sorted BAM
# ============================
echo "removing duplicates for ${prefix}"; date
deduped_bam_dir=${dir}/bam_dedupped
mkdir -p ${deduped_bam_dir}
FINAL_BAM_PREFIX=${prefix}.nodup
FINAL_BAM_FILE="${deduped_bam_dir}/${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${deduped_bam_dir}/${FINAL_BAM_FILE}.bai"
FINAL_BAM_FILE_MAPSTATS="${deduped_bam_dir}/${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}
samtools index ${FINAL_BAM_FILE}
samtools sort -n --threads 10 ${FINAL_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ${FINAL_BAM_FILE_MAPSTATS}

# =============================
# Compute library complexity
# =============================
echo "computing library complexity for ${prefix}"; date
# Sort by name
# convert to bedPE and obtain fragment coordinates
# sort by position and strand
# Obtain unique count statistics
library_complexity_stats_dir=${dir}/library_complexity_stats
mkdir -p ${library_complexity_stats_dir}
PBC_FILE_QC="${library_complexity_stats_dir}/${FINAL_BAM_PREFIX}.pbc.qc"

# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
# col1: total read pairs: total fragments before dedup
# col2: DistinctReadPairs: fragment counts after dedup
# col3: OneReadPair: fragment counts that appeared once
# col4: TwoReadPairs: fragment counts that appeared twice
# col5: NRF = Distinct/total: Non-Redundant Fraction
# col6: PBC1 = OnePair/total: PCR Bottleneck Coefficient 1, >0.9 is great (no bottle necking. 0.5-0.8 is moderate, and 0.8-0.9 is mild.
# col7: PBC2 = OnePair/TwoPair: PCR BottleneckCoefficient2, Ratio between unique and PCR duplicates, >>10 is great
samtools sort -n ${FILT_BAM_FILE} -o ${deduped_bam_dir}/${prefix}.srt.tmp.bam
bedtools bamtobed -bedpe -i ${deduped_bam_dir}/${prefix}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}






























