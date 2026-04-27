#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

module load bedtools
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/heatmap_RNA_round2/TF_peak_overlapping
cd ${dir}

########## fetch QTLs that overlap ChIPseq peaks
cat ${dir}/../reordered_inducedGenes_log2FC1.5_cluster*.txt > gene.txt
enhancer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method4/dictionary_enhancer_updownstreamGenes_all3cts.bed
grep -w -f gene.txt ${enhancer} > ${dir}/existing_enhancers.bed
cat ${dir}/existing_enhancers.bed | cut -f6 | sort | uniq > gene_exist.txt
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
while read gene;do
	grep -w ${gene} ${dir}/existing_enhancers.bed > ${gene}_enhancer.bed
	for f in *.bed.gz;do
		TF=$(echo ${f} | cut -d'_' -f2 | sed 's/-human//g')
		ID=$(echo ${f} | cut -d'_' -f1)
		bedtools intersect -a ${gene}_enhancer.bed -b ${f} | awk -v TF=${TF} -v ID=${ID} '{OFS=FS="\t"}{print $0,TF,ID}' >> enhancer_overlapping_TF_peaks.bed
	done
	echo ${gene}
	rm ${gene}_enhancer.bed
done < ${dir}/gene_exist.txt
mv enhancer_overlapping_TF_peaks.bed ${dir}
cd ${dir}

# clean-up
echo -e "enhancer\tgene\tcelltype\tTF" > enhancer_overlapping_TF_peaks_cleaned.txt
cat enhancer_overlapping_TF_peaks.bed | cut -f4,6,7,9 | grep -v ENCFF585XWV.bed.gz | awk '{OFS=FS="\t"}{print $0,$1"_"$4}' | awk '!seen[$5]++' | awk '{OFS=FS="\t"}{print $1,$2,$3,$4}' >> enhancer_overlapping_TF_peaks_cleaned.txt
#
#Rscript ${dir}/analysis_TF_enhancer_overlap.R # heatmap, venn-diagram, differential enrichment analysis

