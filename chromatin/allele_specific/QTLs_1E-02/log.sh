# written by shuo.shan@umassmed.edu 04/2024

module load bcftools
module load bedtools
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/allele_specific/QTLs_1E-02
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

########## intersect with skin AS SNVs
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENTEX
cat hetSNVs_pooled_AS_DNase.tsv | head -1 > hetSNVs_skin_AS_DNase.txt
cat hetSNVs_pooled_AS_DNase.tsv | grep skin >> hetSNVs_skin_AS_DNase.txt
cat hetSNVs_skin_AS_DNase.txt | awk 'NR>1{OFS=FS="\t"}{if ($17=="1") print $0}' | bedtools sort -i stdin > hetSNVs_skin_AS_DNase_sig.bed
cd ${dir}
ASfile=/pi/manuel.garber-umw/human/skin/eQTLs/literature/ENTEX/hetSNVs_skin_AS_DNase_sig.bed
header1=$(echo -e "#chr\tstart\tend\tID\tREF\tALT")
header2=$(head -1 /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENTEX/hetSNVs_pooled_AS_DNase.tsv)
echo -e "${header1}\t${header2}\tbp_overlap" > QTL_1E-2_overlapping_AlleleSpecific_hetSNVs_DNase_skin.bed
bedtools intersect -wo -a ${dir}/QTL.bed -b ${ASfile} | bedtools sort -i stdin >> QTL_1E-2_overlapping_AlleleSpecific_hetSNVs_DNase_skin.bed
# 29,038

########## intersect with other mark AS hetSNVs
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENTEX
cat hetSNVs_pooled_AS.tsv | awk 'NR>1{if ($17=="1") print $0}' > hetSNVs_pooled_AS_sig.txt
cat hetSNVs_pooled_AS_sig.txt | awk '{print $10}' | sort | uniq > marks.txt
while read mark;do
  cat hetSNVs_pooled_AS_sig.txt | grep -w ${mark} | bedtools sort -i stdin > temp_hetSNVs_pooled_AS_sig.bed
  header1=$(echo -e "#chr\tstart\tend\tID\tREF\tALT")
  header2=$(head -1 /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENTEX/hetSNVs_pooled_AS_DNase.tsv)
  echo -e "${header1}\t${header2}\tbp_overlap" > ${dir}/QTL_1E-2_overlapping_AlleleSpecific_${mark}_hetSNVs_pooledtissue.bed
  bedtools intersect -wo -a ${dir}/QTL.bed -b temp_hetSNVs_pooled_AS_sig.bed |\
	  bedtools sort -i stdin >> ${dir}/QTL_1E-2_overlapping_AlleleSpecific_${mark}_hetSNVs_pooledtissue.bed
  echo "Done with ${mark}"
  rm temp_hetSNVs_pooled_AS_sig.bed
done < marks.txt
cd ${dir}

########## annotate QTLs with all marks if allele specific hetSNVs are detected
# note: output file is QTL_1E-02_annotated_AlleleSpecificMark_summary.bed. only DNase is skin. all other marks are from pooled tissue.
Rscript ${dir}/annotate_QTLs.R

########## clean-up
rm *_*QTL.txt
rm QTL.bed QTL.txt QTL.vcf.gz









######## note#1: this is how to get the reQTLs 
# reQTL cutoff is 1E-05 from rankNormCPM featureSelected model
# fetch MEL reQTL 
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.00001) print $1}' | sort | uniq > MEL_reQTL.txt
# fetch KRT reQTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/KRT_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.00001) print $1}' | sort | uniq > KRT_reQTL.txt
# fetch FRB reQTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/FRB_new/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.00001) print $1}' | sort | uniq > FRB_reQTL.txt
# get unique reQTLs and the vcf file
cat *_reQTL.txt | sort | uniq > reQTL.txt
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
module load bcftools
bcftools view --include ID==@reQTL.txt ${vcf} -Oz -o ${dir}/reQTL.vcf.gz
rm *_reQTL.txt
######## end note#1
