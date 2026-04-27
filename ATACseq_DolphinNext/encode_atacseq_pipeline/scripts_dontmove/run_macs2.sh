#!/bin/bash
# this script in ENCODE ATACseq pipeline performs post-alignment filtering.

# -------------------------------------------------------------------------
# set-up
prefix=$1
dir=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline
deduped_bam_dir=${dir}/bam_dedupped
bedpe_dir=${dir}/BEDPE_TAGALIGN
mkdir -p ${bedpe_dir}
cd ${bedpe_dir}

FINAL_BAM_PREFIX=${prefix}.nodup
FINAL_BAM_FILE=${deduped_bam_dir}/${FINAL_BAM_PREFIX}.bam

source activate ATACseq
module load samtools
module load bedtools
module load openjdk

samtools sort -n ${FINAL_BAM_FILE} -o ${deduped_bam_dir}/${prefix}.srt.tmp.bam
# ================
# Create BEDPE file
# ================
echo "creating BEDPE file for ${prefix}"; date
FINAL_BEDPE_FILE=${bedpe_dir}/${FINAL_BAM_PREFIX}.bedpe.gz
bedtools bamtobed -bedpe -i ${deduped_bam_dir}/${prefix}.srt.tmp.bam | grep -v "chrM" | gzip -nc > ${FINAL_BEDPE_FILE}


# ===================
# Create tagAlign file
# ===================
echo "creating tagAlign file for ${prefix}"; date
FINAL_TA_FILE=${bedpe_dir}/${FINAL_BAM_PREFIX}.tagAlign.gz
zcat ${FINAL_BEDPE_FILE} | \
awk 'BEGIN{OFS="\t"}{
    printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",
    $1,$2,$3,$9,$4,$5,$6,$10
}' | \
awk 'BEGIN{OFS="\t"}{
    if ($6 == "+") {$2 = $2 + 4}
    else if ($6 == "-") {$3 = $3 - 5}
    print $0
}' | gzip -nc > ${FINAL_TA_FILE}


# ===================
# MACS2 peak calling
# ===================
macs2_dir=${dir}/MACS2
mkdir -p ${macs2_dir}
cd ${macs2_dir}
tag=${FINAL_TA_FILE}
gensz="hs"           
pval_thresh=0.01
smooth_window=150
shiftsize=$(( -$smooth_window / 2 ))
chrsz=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
NPEAKS=300000
module load macs2

echo "MACS2 callpeak for ${prefix}"; date
macs2 callpeak \
    -t $tag -f BED -n "$prefix" -g "$gensz" -p $pval_thresh \
   --shift $shiftsize  --extsize $smooth_window --nomodel -B --SPMR --keep-dup all --call-summits

echo "sort narrowPeak for ${prefix}"; date
# Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>
sort -k 8gr,8gr "$prefix"_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | head -n ${NPEAKS} | gzip -nc > ${prefix}.narrowPeak.gz
macs2 bdgcmp -t "$prefix"_treat_pileup.bdg -c "$prefix"_control_lambda.bdg --o-prefix "$prefix" -m FE
# clip peaks to genome bounds
bedClip ${prefix}.narrowPeak.gz $chrsz temp.clipped.bed
sort -k1,1 -k2,2n temp.clipped.bed | gzip -nc > ${prefix}.clipped.narrowPeak.gz






























