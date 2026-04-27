
dir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/test1
cd ${dir}

# get genotype for SNP rs2970686: 5:96916885
f=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/merged.vcf.gz # genotype data before MIS imputation
f=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38.vcf.gz # genotype data post MIS imputation
bcftools view -r chr5:96916885 ${f} | bcftools view -Oz -o snp.vcf.gz
bcftools index snp.vcf.gz

# vcf to bed
bcftools query -f '%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' snp.vcf.gz \
  | awk 'BEGIN {OFS="\t"} { $3 = $3 + 1; print }' > snp.bed
bcftools query -l snp.vcf.gz | paste -sd '\t' - > sample_names.txt
echo -e "CHROM\tSTART\tEND\tID\tREF\tALT\t$(cat sample_names.txt)" > snp.bed.header
cat snp.bed >> snp.bed.header
mv snp.bed.header snp.bed
rm sample_names.txt

genotype_table=snp.vcf.gz
bcftools query -f "%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT{0}[\t%GT]\n" -H ${genotype_table} | head -1 > header
cat header | tr '\t' '\n' | cut -d']' -f2 | tr '\n' '\t'  | sed '$s/\t$/\n/' > header2
cat header2 | sed 's/POS\tPOS/START\tEND/g'  > header3
bcftools query -f "%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT{0}[\t%GT]\n" ${genotype_table} -o temp.filtered.genotype.vcf
cat temp.filtered.genotype.vcf | awk '{OFS="\t"}{print $1,$2,$2+1,$4,$5,$6}' > temp1
cat temp.filtered.genotype.vcf | cut -f7- > temp2
paste temp1 temp2 > temp3
cat temp3 > snp.bed
