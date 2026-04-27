#!/bin/bash
module load bedtools

workingdir=$1
QTLbed=$2
prefix=$3

########## fetch QTLs that overlap ChIPseq peaks
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
rm -f QTL_overlapping_TF_peaks_${prefix}.bed
for f in *.bed.gz;do
        TF=$(echo ${f} | cut -d'_' -f2 | sed 's/-human//g')
        ID=$(echo ${f} | cut -d'_' -f1)
	bedtools intersect -a ${QTLbed} -b ${f} | cut -f1-6 | awk -v TF=${TF} -v ID=${ID} '{OFS=FS="\t"}{print $0,TF,ID}' >> QTL_overlapping_TF_peaks_${prefix}.bed
        echo ${f}
done
mv QTL_overlapping_TF_peaks_${prefix}.bed ${workingdir}
cd ${workingdir}

# print result quickly
wc -l ${workingdir}/QTL_overlapping_TF_peaks_${prefix}.bed

# this create the sorted bed file of high-level summary
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
Rscript /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/annotate_SNP_by_overlapping_TFBSpeaks_compile.R ${workingdir} ${workingdir}/QTL_overlapping_TF_peaks_${prefix}.bed ${prefix} 
