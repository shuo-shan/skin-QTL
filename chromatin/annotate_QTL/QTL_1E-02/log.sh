module load bcftools
module load bedtools
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02
cd ${dir}

########## QTL pval cutoff is 1E-02 from rankNormCPM featureSelected model
# fetch MEL QTL 
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.01) print $1}' | sort | uniq > MEL_reQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($15!="." && $15<0.01) print $1}' | sort | uniq > MEL_PBSeQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($23!="." && $23<0.01) print $1}' | sort | uniq > MEL_IFNeQTL.txt
cat MEL_*QTL.txt | sort | uniq > MEL_QTL.txt # (334,722)
echo -e "QTL\tgene\treQTL_pval\treQTL_beta\treQTL_se\tPBSeQTL_pval\tPBSeQTL_beta\tPBSeQTL_se\tIFNeQTL_pval\tIFNeQTL_beta\tIFNeQTL_se" > MEL_QTL_modeling_result.txt
cat ${f} | grep -w -f MEL_QTL.txt | awk '{OFS=FS="\t"}{print $1,$2,$7,$9,$10,$15,$17,$18,$23,$25,$26}' >> MEL_QTL_modeling_result.txt
# fetch KRT QTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/KRT_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.01) print $1}' | sort | uniq > KRT_reQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($15!="." && $15<0.01) print $1}' | sort | uniq > KRT_PBSeQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($23!="." && $23<0.01) print $1}' | sort | uniq > KRT_IFNeQTL.txt
cat KRT_*QTL.txt | sort | uniq > KRT_QTL.txt # (320,971)
echo -e "QTL\tgene\treQTL_pval\treQTL_beta\treQTL_se\tPBSeQTL_pval\tPBSeQTL_beta\tPBSeQTL_se\tIFNeQTL_pval\tIFNeQTL_beta\tIFNeQTL_se" > KRT_QTL_modeling_result.txt
cat ${f} | grep -w -f KRT_QTL.txt | awk '{OFS=FS="\t"}{print $1,$2,$7,$9,$10,$15,$17,$18,$23,$25,$26}'  >> KRT_QTL_modeling_result.txt
# fetch FRB QTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/FRB_new/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.01) print $1}' | sort | uniq > FRB_reQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($15!="." && $15<0.01) print $1}' | sort | uniq > FRB_PBSeQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($23!="." && $23<0.01) print $1}' | sort | uniq > FRB_IFNeQTL.txt
cat FRB_*QTL.txt | sort | uniq > FRB_QTL.txt # (323,391)
echo -e "QTL\tgene\treQTL_pval\treQTL_beta\treQTL_se\tPBSeQTL_pval\tPBSeQTL_beta\tPBSeQTL_se\tIFNeQTL_pval\tIFNeQTL_beta\tIFNeQTL_se" > FRB_QTL_modeling_result.txt
cat ${f} | grep -w -f FRB_QTL.txt | awk '{OFS=FS="\t"}{print $1,$2,$7,$9,$10,$15,$17,$18,$23,$25,$26}' >> FRB_QTL_modeling_result.txt
# get unique QTLs and the vcf file and bed file
cat MEL_QTL.txt KRT_QTL.txt FRB_QTL.txt | sort | uniq > QTL.txt #(691,034)
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
bcftools view -m2 -M2 -v snps --include ID==@QTL.txt ${vcf} -Oz -o ${dir}/QTL.vcf.gz
bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\t%ID\n' ${dir}/QTL.vcf.gz | awk '{OFS=FS="\t"}{print $1,$2-1,$2,$6,$4,$5}' | bedtools sort -i stdin > ${dir}/QTL.bed

# allele-specific
f1=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/allele_specific/QTLs_1E-02/QTL_1E-02_annotated_AlleleSpecificMark.bed
# ANNOVAR
f2=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/ANNOVAR/QTLs_1E-02/QTL_1E-02_annotated_ANNOVAR.bed
# TF_peak
f3=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/TF_peak_overlapping/QTL_1E-02/QTL_overlapping_TF_peaks_collapsed.bed
# TF motif
f4=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/TF_motif_overlapping/QTL_1E-02/QTL_1E-02_overlapping_TFmotif.bed
# cRE and dynamics
f5=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/promoter_enhancer_overlapping/QTLs_1E-02/QTL_overlapping_CRE.bed
# beta-comparison modeling results

# tables are joined by Rscript join_tables.R 

# results are explored in Rscript explore_results.R



# compiled tables are saved as .txt files.
# significant cutoff is: either p.linearModel < 1E-05 or p.permute < 0.001
for celltype in MEL KRT FRB;do
	echo "starting with celltype: ${celltype}";date
	f=compiled_table_${celltype}_reQTL.txt
	cat ${f} | cut -f2 | awk 'NR>1' | sort | uniq > this_genes.txt
	while read gene;do
		mkdir -p /pi/manuel.garber-umw/human/skin/eQTLs/website/data/plots/${gene}
		cat ${f} | awk -v g=${gene} 'NR==1{print $0}NR>1{if ($2==g) print $0}' > /pi/manuel.garber-umw/human/skin/eQTLs/website/data/plots/${gene}/compiled_table_${celltype}_reQTL.txt
	done < this_genes.txt
	rm this_genes.txt
	echo "done with celltype: ${celltype}";date
done
