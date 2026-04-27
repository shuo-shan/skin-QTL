#!/bin/bash
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=32000]"
#BSUB -W 01:00
#BSUB -q short
#BSUB -J mergeCov
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/log/mergeCovAcrossSamples_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/log/mergeCovAcrossSamples_%J.err"

# ------------------
# merge all sample-level covarege of masterpeaks
# ------------------
BAM_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped
COL_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/cols
MASTER=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/master_peaks.bed
OUT=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/master_peaks_multicov.txt

BAM_LIST=($(ls ${BAM_DIR}/*.nodup.bam | sort))

# Build header
HEADER="chr\tstart\tend\tpeak_name"
for bam in "${BAM_LIST[@]}"; do
    HEADER="${HEADER}\t$(basename ${bam} .nodup.bam)"
done

# Collect count columns in the same sorted order
COL_FILES=()
for bam in "${BAM_LIST[@]}"; do
    SAMPLE=$(basename ${bam} .nodup.bam)
    COL_FILES+=("${COL_DIR}/${SAMPLE}.counts.txt")
done

# Sanity check: all files present?
for f in "${COL_FILES[@]}"; do
    if [ ! -f "$f" ]; then echo "MISSING: $f"; exit 1; fi
done

echo -e "${HEADER}" > ${OUT}
paste ${MASTER} "${COL_FILES[@]}" >> ${OUT}

echo "Done: $(wc -l ${OUT} | awk '{print $1}') lines written"; date

# ------------------
# create bed file
# ------------------
module load htslib/1.16
cd /pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks
tail -n +2 master_peaks_multicov.txt > master_peaks_multicov.bed
bgzip master_peaks_multicov.bed
tabix -p bed master_peaks_multicov.bed.gz

# ------------------
# QC metric
# ------------------
MASTER=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/master_peaks.bed
MULTICOV=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/master_peaks_multicov.txt
OUT_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/QC
mkdir -p ${OUT_DIR}
echo "=== Peak size distribution ==="
awk 'BEGIN{OFS="\t"} {len=$3-$2; print $0, len}' ${MASTER} \
    | sort -k5,5rn \
    > ${OUT_DIR}/peaks_with_size.bed

# Largest peak
echo "Largest peak:"
head -1 ${OUT_DIR}/peaks_with_size.bed

# Summary stats (min, max, mean, median)
awk '{print $3-$2}' ${MASTER} | sort -n | awk '
    BEGIN{min=999999; max=0; sum=0; n=0}
    {
        a[n]=$1; sum+=$1; n++;
        if($1<min) min=$1;
        if($1>max) max=$1;
    }
    END{
        mean=sum/n;
        median=(n%2==0) ? (a[int(n/2)-1]+a[int(n/2)])/2 : a[int(n/2)];
        printf "N peaks:\t%d\nMin size:\t%d bp\nMax size:\t%d bp\nMean size:\t%.1f bp\nMedian size:\t%.1f bp\nTotal bases:\t%d bp\n", n, min, max, mean, median, sum
    }' | tee ${OUT_DIR}/peak_size_summary.txt


# Size distribution histogram buckets
echo ""
echo "=== Size distribution (bp) ==="
awk '{print $3-$2}' ${MASTER} | awk '
    {
        if($1<=300) b["<=300"]++;
        else if($1<=500) b["301-500"]++;
        else if($1<=1000) b["501-1000"]++;
        else if($1<=2000) b["1001-2000"]++;
        else b[">2000"]++;
    }
    END{
        print "<=300\t" b["<=300"]+0;
        print "301-500\t" b["301-500"]+0;
        print "501-1000\t" b["501-1000"]+0;
        print "1001-2000\t" b["1001-2000"]+0;
        print ">2000\t" b[">2000"]+0;
    }' | tee ${OUT_DIR}/size_distribution.txt


# Per-chromosome peak counts
echo ""
echo "=== Peaks per chromosome ==="
awk '{print $1}' ${MASTER} | sort | uniq -c | sort -k2,2V \
    | awk '{print $2"\t"$1}' | tee ${OUT_DIR}/peaks_per_chrom.txt

# Coverage completeness: how many peaks have zero counts across ALL samples?
echo ""
echo "=== Zero-count peaks (likely noise) ==="
awk 'NR>1 {
    zero=1;
    for(i=5;i<=NF;i++) if($i>0){zero=0; break}
    if(zero) zero_count++
} END{print zero_count " peaks with zero counts in all samples"}' ${MULTICOV} \
    | tee -a ${OUT_DIR}/peak_size_summary.txt





































