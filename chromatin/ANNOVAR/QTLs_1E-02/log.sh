# written by shuo.shan@umassmed.edu 04/2024

module load bcftools
module load bedtools
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/ANNOVAR/QTLs_1E-02
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

# compile to ANNOVAR format
bcftools query -f '%CHROM %POS %POS %REF %ALT %ID\n' ${dir}/QTL.vcf.gz > ${dir}/QTL_annovar_input.txt

# run ANNOVAR
# output file: ${dir}/reQTL.hg38_multianno.txt. 
# header annotation is ${dir}/annovar_output_header_notes.txt
annovarDir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/ANNOVAR/annovar
perl ${annovarDir}/table_annovar.pl ${dir}/QTL_annovar_input.txt ${annovarDir}/humandb/ -buildver hg38 -out QTL -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -polish -xref ${annovarDir}/example/gene_xref.txt

# classify based on Func.refGene
cat ${dir}/QTL.hg38_multianno.txt | awk 'NR>1' | cut -f6 | sort | uniq -c > summary_Func.refGene.txt  
cat ${dir}/QTL.hg38_multianno.txt | awk 'NR>1' | cut -f9 | sort | uniq -c > summary_ExonicFunc.refGene.txt


# investigate exonic variants
cat ${dir}/QTL.hg38_multianno.txt | grep exonic | grep nonsynonymous | cut -f21 | sort | uniq > temp_variants.txt
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.00001) print $0}' | grep -w -f temp_variants.txt | cut -f1,2 > temp_MEL.txt
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/KRT_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.00001) print $0}' | grep -w -f temp_variants.txt | cut -f1,2 > temp_KRT.txt
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/FRB_new/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.00001) print $0}' | grep -w -f temp_variants.txt | cut -f1,2 > temp_FRB.txt
cat temp_MEL.txt temp_KRT.txt temp_FRB.txt | awk '{print $1"_"$2}' > temp_snp_gene_pairs.txt # use R to make plots
cat ${dir}/QTL.hg38_multianno.txt | grep exonic | grep nonsynonymous | cut -f1,2,4,5,6,7,9,11,21


# compile a summary file on refGene
#echo -e "chr\tstart\tend\tID\tREF\tALT\tANNOVAR_Func_Gene_ExonicFunc_refGene" > QTL_1E-02_annotated_ANNOVAR.bed
#cat ${dir}/QTL.hg38_multianno.txt | awk '{OFS=FS="\t"}NR>1{print $1,$2-1,$3,$21,$4,$5,$6"_"$7"_"$9}' | sed 's/ //g' >> QTL_1E-02_annotated_ANNOVAR.bed




# double check the SNP information agree
cat QTL_1E-02_annotated_ANNOVAR.bed | awk 'NR>1' | cut -f2 > temp1.txt
cat ${dir}/QTL.bed | cut -f2 > temp2.txt
comm -1 -2 <(sort temp1.txt) <(sort temp2.txt) | grep -v '^$' | wc -l
