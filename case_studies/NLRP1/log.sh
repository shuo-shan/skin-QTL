#!/bin/bash

# set-up
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/case_studies/NLRP1
g=NLRP1

vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz
vcf_subset=NLRP1.hap.tags.vcf.gz
samples=samples.txt
header=header
body=body
geno=temp.snp.txt
geno_num=NLRP1.hap.tags.txt

cd ${DIR}

# ---------- get gene, SNP, metadata, and modeling stats info ----------- #
# get gene CPM
cpm_all_f=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt
head -1 ${cpm_all_f} > ${DIR}/cpm.txt
awk -v gene=${g} '{OFS=FS="\t"}{if ($NF==gene) print $0}' ${cpm_all_f} >> ${DIR}/cpm.txt
echo "got gene CPM"

# Extract Haplotype 1 and 2A SNPs from our VCF
echo "==== Subsetting genotype VCF for ${snp} ===="; date
module load bcftools

# 1) extract haplotype defining SNPs from VCF file
cut -f1 NLRP1_haplotype_definitions.tsv | tail -n +2 > NLRP1_hap_rsids.txt
bcftools view -i 'ID=@NLRP1_hap_rsids.txt' \
  -Oz -o ${vcf_subset} ${vcf}

# 2) clean sample names
bcftools query -l "$vcf_subset" | cut -d'_' -f1 | sed 's/F0/F/g' | sed 's/skineQTL-//g' > "$samples"

# 3) build header
printf "CHROM\tPOS\tID\tREF\tALT\t" > "$header"
paste -sd '\t' "$samples" >> "$header"

# 4) extract GT matrix
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$vcf_subset" > "$body"
cat "$header" "$body" > "$geno"

# 5) convert genotypes to numeric
sed -E 's/0\/0/0/g; s/0\/1/1/g; s/1\/0/1/g; s/1\/1/2/g; s/\.\/\./NA/g' "$geno" > "$geno_num"

# 6) cleanup intermediate files (leave only numeric genotype)
rm "$vcf_subset" "$samples" "$header" "$body" "$geno"

echo "got SNP genotype"

