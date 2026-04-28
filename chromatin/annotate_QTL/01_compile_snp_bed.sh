#!/bin/bash
# read: a SNP ID
# output: a SNP bed file

workingdir=$1
SNPID=$2

module load bcftools
module load bedtools

queryBedF=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/SNPs_near_TSS.bed
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz

cd ${workingdir}
echo ${SNPID} > temp_${SNPID}.txt
bcftools view --include ID==@temp_${SNPID}.txt ${vcf} -Oz -o ${workingdir}/${SNPID}.vcf.gz
bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\t%ID\n' ${workingdir}/${SNPID}.vcf.gz | awk '{OFS=FS="\t"}{print $1,$2-1,$2,$6,$4,$5}' | bedtools sort -i stdin > ${workingdir}/${SNPID}.bed

rm -f ${workingdir}/temp_${SNPID}.txt
rm -f ${workingdir}/${SNPID}.vcf.gz

echo "wrote to ${workingdir}/${SNPID}.bed"; date

