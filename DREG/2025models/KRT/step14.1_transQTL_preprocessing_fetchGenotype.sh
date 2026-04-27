#!/bin/bash

module load bcftools

ct=$1
cond=$2
QTLtype=$3
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}

# ---------------------------------- #
# compile genotype matrix for QTLs
# ---------------------------------- #
echo "==== Subsetting genotype VCF for ${ct}_${cond}_${QTLtype}_SNPs ===="; date
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz
vcf_subset=${DIR}/transQTL/QTL_tags/${ct}_${cond}_${QTLtype}_SNPs.vcf.gz
samples=${DIR}/transQTL/QTL_tags/${ct}_${cond}_${QTLtype}_samples.txt
header=${DIR}/transQTL/QTL_tags/${ct}_${cond}_${QTLtype}_header
body=${DIR}/transQTL/QTL_tags/${ct}_${cond}_${QTLtype}_body
geno=${DIR}/transQTL/QTL_tags/${ct}_${cond}_${QTLtype}_temp.snp.txt
geno_num=${DIR}/transQTL/QTL_tags/${ct}_${cond}_${QTLtype}_SNPs_genotype.txt

# 1) collect SNP IDs from this chunk
QTLtags=${DIR}/transQTL/QTL_tags/${ct}_${cond}_${QTLtype}_QTLtags.txt
SNP_list=${DIR}/transQTL/QTL_tags/${ct}_${cond}_${QTLtype}_SNPs.txt
cat ${QTLtags} | awk 'NR>1' | cut -f1 > ${SNP_list}

# 2) extract VCF records
bcftools view --include ID==@${SNP_list} ${vcf} -Oz -o ${vcf_subset}

# 3) clean sample names
bcftools query -l "$vcf_subset" | cut -d'_' -f1 | sed 's/F0/F/g' | sed 's/skineQTL-//g' > "$samples"

# 4) build header
printf "CHROM\tPOS\tID\tREF\tALT\t" > "$header"
paste -sd '\t' "$samples" >> "$header"

# 5) extract GT matrix
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$vcf_subset" > "$body"
cat "$header" "$body" > "$geno"

# 6) convert genotypes to numeric
sed -E 's/0\/0/0/g; s/0\/1/1/g; s/1\/0/1/g; s/1\/1/2/g; s/\.\/\./NA/g' "$geno" > "$geno_num"

# 7) cleanup intermediate files (leave only numeric genotype)
rm "${SNP_list}" "$vcf_subset" "$samples" "$header" "$body" "$geno"

echo "got SNP genotype"
