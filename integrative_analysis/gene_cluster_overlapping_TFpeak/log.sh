#!/bin/bash
# written by shuo.shan@umassmed.edu 05/2024
# overlap commonly induced genes' promoters with ENCODE skin tissue TF ChIPseq sig peaks

dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/gene_cluster_overlapping_TFpeak

module load bedtools
cd ${dir}

### commonly induced promoters
# intersect with ENCODE TF peaks
promoterF=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/XSTREME/commonly_induced_genes_promoter/commonly_induced_promoters.bed
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
rm promoter_overlapping_TF_peaks.bed
for f in *.bed.gz;do
        TF=$(echo ${f} | cut -d'_' -f2 | sed 's/-human//g')
        ID=$(echo ${f} | cut -d'_' -f1)
        bedtools intersect -a ${promoterF} -b ${f} | cut -f1-6 | awk -v TF=${TF} -v ID=${ID} '{OFS=FS="\t"}{print $0,TF,ID}' >> promoter_overlapping_TF_peaks.bed
        echo ${f}
done
mv promoter_overlapping_TF_peaks.bed ${dir}/commonly_induced_genes_promoter_overlapping_TF_peaks.bed
cd ${dir}

### commonly induced promoters
# intersect with ENCODE TF peaks
promoterF=

# check which genes are regulated by which TF
# overlap_with_TF_peak.R
# output is: ENCODE_TF_peak_intersecting_promoters_of_commonly_induced_genes_summary.txt

