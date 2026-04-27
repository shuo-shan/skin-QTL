# written by shuo.shan@umassmed.edu 04/2024
# for each cell type, ask this question for each QTL: 
# 1. CRE_overlap: does it overlap a promoter/enhancer (promoter/enhancer)? 
# 3. CRE_dynamic: what is the dynamic of the promoter/enhancer it overlaps/near? (open/close, gain/loseK27ac)
# promoter: < 500bp of expressed gene TSS
# enhancer: < 300Kbp of expressed gene TSS
# dynamic: padj < 0.1, log2FC > or < 0
module load bcftools
module load bedtools
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/promoter_enhancer_overlapping/QTLs_1E-02
cd ${dir}

########## QTL pval cutoff is 1E-02 from rankNormCPM featureSelected model
# fetch MEL QTL 
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.01) print $1}' | sort | uniq > MEL_reQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($15!="." && $15<0.01) print $1}' | sort | uniq > MEL_PBSeQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($23!="." && $23<0.01) print $1}' | sort | uniq > MEL_IFNeQTL.txt
cat MEL_*QTL.txt | sort | uniq > MEL_QTL.txt # (334,722)
# fetch KRT QTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/KRT_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.01) print $1}' | sort | uniq > KRT_reQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($15!="." && $15<0.01) print $1}' | sort | uniq > KRT_PBSeQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($23!="." && $23<0.01) print $1}' | sort | uniq > KRT_IFNeQTL.txt
cat KRT_*QTL.txt | sort | uniq > KRT_QTL.txt # (320,971)
# fetch FRB QTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/FRB_new/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.01) print $1}' | sort | uniq > FRB_reQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($15!="." && $15<0.01) print $1}' | sort | uniq > FRB_PBSeQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($23!="." && $23<0.01) print $1}' | sort | uniq > FRB_IFNeQTL.txt
cat FRB_*QTL.txt | sort | uniq > FRB_QTL.txt # (323,391)
# get unique QTLs and the vcf file and bed file
cat MEL_QTL.txt KRT_QTL.txt FRB_QTL.txt | sort | uniq > QTL.txt #(691,034)
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
bcftools view --include ID==@QTL.txt ${vcf} -Oz -o ${dir}/QTL.vcf.gz
bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\t%ID\n' ${dir}/QTL.vcf.gz | awk '{OFS=FS="\t"}{print $1,$2-1,$2,$6,$4,$5}' | bedtools sort -i stdin > ${dir}/QTL.bed

########## get promoters, enhancers, and dynamics on these regions for each celltype
########## FRB
address=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/merged_coverage_ATACseq_allcts_allregions_1kbp.bed
ATACdynamic=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_FRB.txt
ChIPdynamic=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DESeq2/DESeq2_results_H3K27acChIPseq_peaks_1kb_FRB.txt

promoter=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_gene_promoter_links_FRB.bed
enhancer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_closest_gene_enhancer_links_FRB.bed

cat ${address} | awk 'NR>1' | cut -f1-4 | sort -k4,4V > CRE_1kbp_peak.bed

cat ${promoter} | awk 'NR>1' | cut -f5 | sort | uniq | awk '{OFS="\t"}{print $1,"promoter"}' > temp_promoter.txt
cat ${enhancer} | awk 'NR>1' | cut -f5 | sort | uniq | awk '{OFS="\t"}{print $1,"enhancer"}' > temp_enhancer.txt
cat temp_promoter.txt temp_enhancer.txt | cut -f1 | sort | uniq > peak_to_keep.txt
cat CRE_1kbp_peak.bed |  grep -v -w -f peak_to_keep.txt  | cut -f4 | awk '{OFS="\t"}{print $1,"nonCREpeak"}' > temp_na.txt
cat temp_promoter.txt temp_enhancer.txt temp_na.txt | sort -k1,1V > CRE_1kbp_info.txt

cat ${ATACdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3>0) print $1,"open"}' > temp1.txt
cat ${ATACdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3<0) print $1,"close"}' > temp2.txt
cat temp1.txt temp2.txt | cut -f1 | sort | uniq > peak_to_keep.txt
cat ${ATACdynamic} | grep -v -w -f peak_to_keep.txt | awk 'NR>1' | awk '{OFS="\t"}{print $1,"NA"}' > temp3.txt
cat temp1.txt temp2.txt temp3.txt | sort -k1,1V > CRE_1kbp_dynamic_atac.txt

cat ${ChIPdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3>0) print $1,"gainK27ac"}' > temp1.txt
cat ${ChIPdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3<0) print $1,"loseK27ac"}' > temp2.txt
cat temp1.txt temp2.txt | cut -f1 | sort | uniq > peak_to_keep.txt
cat ${ChIPdynamic} | grep -v -w -f peak_to_keep.txt | awk 'NR>1' | awk '{OFS="\t"}{print $1,"NA"}' > temp3.txt
cat temp1.txt temp2.txt temp3.txt | sort -k1,1V > CRE_1kbp_dynamic_chip.txt

join -1 4 -2 1 CRE_1kbp_peak.bed CRE_1kbp_info.txt | sort -k1,1V > temp1.txt
join -1 1 -2 1 CRE_1kbp_dynamic_atac.txt CRE_1kbp_dynamic_chip.txt | awk '{OFS="\t"}{print $1,$2"_"$3}' | sort -k1,1V > temp_joined.txt
join -1 1 -2 1 temp1.txt temp_joined.txt | awk '{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' | sort -k4,4V > CRE_dynamic_FRB.bed
rm temp*.txt peak_to_keep.txt CRE_1kbp_peak.bed CRE_1kbp_info.txt CRE_1kbp_dynamic_atac.txt CRE_1kbp_dynamic_chip.txt 

########## KRT
address=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/merged_coverage_ATACseq_allcts_allregions_1kbp.bed
ATACdynamic=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_KRT.txt
ChIPdynamic=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DESeq2/DESeq2_results_H3K27acChIPseq_peaks_1kb_KRT.txt
promoter=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_gene_promoter_links_KRT.bed
enhancer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_closest_gene_enhancer_links_KRT.bed

cat ${address} | awk 'NR>1' | cut -f1-4 | sort -k4,4V > CRE_1kbp_peak.bed

cat ${promoter} | awk 'NR>1' | cut -f5 | sort | uniq | awk '{OFS="\t"}{print $1,"promoter"}' > temp_promoter.txt
cat ${enhancer} | awk 'NR>1' | cut -f5 | sort | uniq | awk '{OFS="\t"}{print $1,"enhancer"}' > temp_enhancer.txt
cat temp_promoter.txt temp_enhancer.txt | cut -f1 | sort | uniq > peak_to_keep.txt
cat CRE_1kbp_peak.bed |  grep -v -w -f peak_to_keep.txt  | cut -f4 | awk '{OFS="\t"}{print $1,"nonCREpeak"}' > temp_na.txt
cat temp_promoter.txt temp_enhancer.txt temp_na.txt | sort -k1,1V > CRE_1kbp_info.txt

cat ${ATACdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3>0) print $1,"open"}' > temp1.txt
cat ${ATACdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3<0) print $1,"close"}' > temp2.txt
cat temp1.txt temp2.txt | cut -f1 | sort | uniq > peak_to_keep.txt
cat ${ATACdynamic} | grep -v -w -f peak_to_keep.txt | awk 'NR>1' | awk '{OFS="\t"}{print $1,"NA"}' > temp3.txt
cat temp1.txt temp2.txt temp3.txt | sort -k1,1V > CRE_1kbp_dynamic_atac.txt

cat ${ChIPdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3>0) print $1,"gainK27ac"}' > temp1.txt
cat ${ChIPdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3<0) print $1,"loseK27ac"}' > temp2.txt
cat temp1.txt temp2.txt | cut -f1 | sort | uniq > peak_to_keep.txt
cat ${ChIPdynamic} | grep -v -w -f peak_to_keep.txt | awk 'NR>1' | awk '{OFS="\t"}{print $1,"NA"}' > temp3.txt
cat temp1.txt temp2.txt temp3.txt | sort -k1,1V > CRE_1kbp_dynamic_chip.txt

join -1 4 -2 1 CRE_1kbp_peak.bed CRE_1kbp_info.txt | sort -k1,1V > temp1.txt
join -1 1 -2 1 CRE_1kbp_dynamic_atac.txt CRE_1kbp_dynamic_chip.txt | awk '{OFS="\t"}{print $1,$2"_"$3}' | sort -k1,1V > temp_joined.txt
join -1 1 -2 1 temp1.txt temp_joined.txt | awk '{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' | sort -k4,4V > CRE_dynamic_KRT.bed
rm temp*.txt peak_to_keep.txt CRE_1kbp_peak.bed CRE_1kbp_info.txt CRE_1kbp_dynamic_atac.txt CRE_1kbp_dynamic_chip.txt

########## MEL
address=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/merged_coverage_ATACseq_allcts_allregions_1kbp.bed
ATACdynamic=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_MEL.txt
ChIPdynamic=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DESeq2/DESeq2_results_H3K27acChIPseq_peaks_1kb_MEL.txt
promoter=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_gene_promoter_links_MEL.bed
enhancer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_closest_gene_enhancer_links_MEL.bed

cat ${address} | awk 'NR>1' | cut -f1-4 | sort -k4,4V > CRE_1kbp_peak.bed

cat ${promoter} | awk 'NR>1' | cut -f5 | sort | uniq | awk '{OFS="\t"}{print $1,"promoter"}' > temp_promoter.txt
cat ${enhancer} | awk 'NR>1' | cut -f5 | sort | uniq | awk '{OFS="\t"}{print $1,"enhancer"}' > temp_enhancer.txt
cat temp_promoter.txt temp_enhancer.txt | cut -f1 | sort | uniq > peak_to_keep.txt
cat CRE_1kbp_peak.bed |  grep -v -w -f peak_to_keep.txt  | cut -f4 | awk '{OFS="\t"}{print $1,"nonCREpeak"}' > temp_na.txt
cat temp_promoter.txt temp_enhancer.txt temp_na.txt | sort -k1,1V > CRE_1kbp_info.txt

cat ${ATACdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3>0) print $1,"open"}' > temp1.txt
cat ${ATACdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3<0) print $1,"close"}' > temp2.txt
cat temp1.txt temp2.txt | cut -f1 | sort | uniq > peak_to_keep.txt
cat ${ATACdynamic} | grep -v -w -f peak_to_keep.txt | awk 'NR>1' | awk '{OFS="\t"}{print $1,"NA"}' > temp3.txt
cat temp1.txt temp2.txt temp3.txt | sort -k1,1V > CRE_1kbp_dynamic_atac.txt

cat ${ChIPdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3>0) print $1,"gainK27ac"}' > temp1.txt
cat ${ChIPdynamic} | awk 'NR>1' | awk '$7!=""' | awk '{OFS="\t"}{if ($7<0.1 && $3<0) print $1,"loseK27ac"}' > temp2.txt
cat temp1.txt temp2.txt | cut -f1 | sort | uniq > peak_to_keep.txt
cat ${ChIPdynamic} | grep -v -w -f peak_to_keep.txt | awk 'NR>1' | awk '{OFS="\t"}{print $1,"NA"}' > temp3.txt
cat temp1.txt temp2.txt temp3.txt | sort -k1,1V > CRE_1kbp_dynamic_chip.txt

join -1 4 -2 1 CRE_1kbp_peak.bed CRE_1kbp_info.txt | sort -k1,1V > temp1.txt
join -1 1 -2 1 CRE_1kbp_dynamic_atac.txt CRE_1kbp_dynamic_chip.txt | awk '{OFS="\t"}{print $1,$2"_"$3}' | sort -k1,1V > temp_joined.txt
join -1 1 -2 1 temp1.txt temp_joined.txt | awk '{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' | sort -k4,4V > CRE_dynamic_MEL.bed
rm temp*.txt peak_to_keep.txt CRE_1kbp_peak.bed CRE_1kbp_info.txt CRE_1kbp_dynamic_atac.txt CRE_1kbp_dynamic_chip.txt


##################### overlap QTL to CREs
join -1 4 -2 4 CRE_dynamic_FRB.bed CRE_dynamic_KRT.bed | sort -k1,1V | awk '{OFS="\t"}{print $2,$3,$4,$1,$5,$6,$10,$11}' > temp1.bed
join -1 4 -2 4 temp1.bed CRE_dynamic_MEL.bed | sort -k1,1V | awk '{OFS="\t"}{print $2,$3,$4,$1,$5,$6,$7,$8,$12,$13}' > CRE_dynamics_3cts.bed

mkdir log
jobname=join
split QTL.bed -l 2000 splitQTL_
for f in splitQTL_*;do
  echo -e "bash script.sh ${f}" >> commands.txt
done

while read c;do
        echo ${c} | bsub -W 04:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=800] -q long -e "${dir}/log/join_%J%I.err" -o "${dir}/log/join_%J%I.out"
done < ${dir}/commands.txt

echo -e "chr\tstart\tend\tID\tREF\tALT\tcRE_chr\tcRE_start\tcRE_end\tcRE_name\tcRE_overlap_FRB\tcRE_dynamic_FRB\tcRE_overlap_KRT\tcRE_dynamic_KRT\tcRE_overlap_MEL\tcRE_dynamic_MEL" > header
cat QTL_overlapping_CRE.txt >> header

mv header QTL_overlapping_CRE.bed
