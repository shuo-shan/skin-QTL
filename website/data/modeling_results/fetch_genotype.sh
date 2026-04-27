#!/bin/bash
# written by Crystal Shan Nov 2023

chr=$1
pos=$2
vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz

module load bcftools
echo ${chr}$'\t'${pos} > fetch_genotype_temp1.txt

cat /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/genotype_file_header.txt | tr '\t' '\n' > fetch_genotype_temp2.txt
bcftools view -O v -R fetch_genotype_temp1.txt ${vcfF} | bcftools query -f '%CHROM %POS %ID %REF %ALT[ %GT]\n' | tr ' ' '\n' | awk 'NF' > fetch_genotype_temp3.txt

paste fetch_genotype_temp2.txt fetch_genotype_temp3.txt
rm fetch_genotype_temp*.txt
