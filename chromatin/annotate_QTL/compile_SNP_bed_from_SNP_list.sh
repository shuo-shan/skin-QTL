#!/bin/bash
# read: list of SNPs, one per row
# output: SNP bed file

workingdir=$1
SNPlist=$2
prefix=$3

module load bcftools
module load bedtools

queryBedF=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/SNPs_near_TSS.bed
cd ${workingdir}
cp ${SNPlist} ./temp_SNPlist_${prefix}.txt

vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz
bcftools view --include ID==@temp_SNPlist_${prefix}.txt ${vcf} -Oz -o ${workingdir}/QTL_${prefix}.vcf.gz
bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\t%ID\n' ${workingdir}/QTL_${prefix}.vcf.gz | awk '{OFS=FS="\t"}{print $1,$2-1,$2,$6,$4,$5}' | bedtools sort -i stdin > ${workingdir}/QTL_${prefix}.bed

rm -f ${workingdir}/temp_SNPlist_${prefix}.txt
rm -f ${workingdir}/QTL_${prefix}.vcf.gz

echo "wrote to ${workingdir}/QTL_${prefix}.bed"; date

