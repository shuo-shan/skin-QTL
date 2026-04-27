module load bedtools
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/heatmap_RNA_round2/TF_peak_overlapping
cd ${dir}

########## fetch QTLs that overlap ChIPseq peaks
cat ${dir}/../reordered_inducedGenes_log2FC1.5_cluster*.txt > gene.txt
promoter=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_gene_promoter_links_allcts.bed
grep -w -f gene.txt ${promoter} > ${dir}/existing_promoters.bed
cat ${dir}/existing_promoters.bed | cut -f6 | sort | uniq > gene_exist.txt
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
while read gene;do
	grep -w ${gene} ${dir}/existing_promoters.bed > ${gene}_promoter.bed
	for f in *.bed.gz;do
		TF=$(echo ${f} | cut -d'_' -f2 | sed 's/-human//g')
		ID=$(echo ${f} | cut -d'_' -f1)
		bedtools intersect -a ${gene}_promoter.bed -b ${f} | awk -v TF=${TF} -v ID=${ID} '{OFS=FS="\t"}{print $0,TF,ID}' >> promoter_overlapping_TF_peaks.bed
	done
	echo ${gene}
	rm ${gene}_promoter.bed
done < ${dir}/gene_exist.txt
mv promoter_overlapping_TF_peaks.bed ${dir}
cd ${dir}

# clean-up
echo -e "promoter\tgene\tcelltype\tATACdynamic\tH3K27acdynamic\tTF" > promoter_overlapping_TF_peaks_cleaned.txt
cat promoter_overlapping_TF_peaks.bed | cut -f4,6,8,12,16,17 | grep -v ENCFF585XWV.bed.gz | awk '{OFS=FS="\t"}{print $0,$1"_"$6}' | awk '!seen[$7]++' | awk '{OFS=FS="\t"}{print $1,$2,$3,$4,$5,$6}' >> promoter_overlapping_TF_peaks_cleaned.txt

Rscript ${dir}/analysis_TF_promoter_overlap.R # heatmap, venn-diagram, differential enrichment analysis

