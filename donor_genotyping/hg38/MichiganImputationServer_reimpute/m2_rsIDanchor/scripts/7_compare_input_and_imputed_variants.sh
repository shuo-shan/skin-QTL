#!/bin/bash
#BSUB -J combine_imputed_vcfs[1-23]
#BSUB -R "rusage[mem=120000]"
#BSUB -o step7_compare_pre_vs_post_MIS_%I.out
#BSUB -e step7_compare_pre_vs_post_MIS_%I.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1

# Load necessary modules
module load bcftools
module load htslib

DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor
CHROM_LIST=({1..22} X)
chr_name="chr${CHROM_LIST[$((LSB_JOBINDEX-1))]}"
echo "Processing ${chr_name}..."

# -------------------------------------------------------------
# -------- get original preMIS vcf for this chromosome --------
# -------------------------------------------------------------
# subset hg38 vcf to this chromosome
VCF_HG38_PREFIX=preMIS_hg38
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/merged_lowGPreset.vcf.gz
bcftools view -r ${chr_name} ${vcf} -Oz -o ${DIR}/data/${VCF_HG38_PREFIX}_${chr_name}.vcf.gz
bcftools view -H ${DIR}/data/${VCF_HG38_PREFIX}_${chr_name}.vcf.gz | cut -f3 | grep "rs" | tr ';' '\n' | sort | uniq > ${DIR}/data/${VCF_HG38_PREFIX}_${chr_name}_variants.txt

# -------------------------------------------------------------
# -------- get final postMIS vcf for this chr -----------------
# -------------------------------------------------------------
postMIS_vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/${chr_name}_imputed_combined_hg38_rsID.vcf.gz
#bcftools view -H ${postMIS_vcf} | grep "rs" | cut -f3 | sort | uniq > ${DIR}/data/postMIS_hg38_${chr_name}_variants.txt


# Comparison to get new variants from imputation
preMIS_only=$(comm -23 ${DIR}/data/${VCF_HG38_PREFIX}_${chr_name}_variants.txt ${DIR}/data/postMIS_hg38_${chr_name}_variants.txt | wc -l | cut -d' ' -f1)
common=$(comm -12 ${DIR}/data/${VCF_HG38_PREFIX}_${chr_name}_variants.txt ${DIR}/data/postMIS_hg38_${chr_name}_variants.txt | wc -l | cut -d' ' -f1)
postMIS_only=$(comm -13 ${DIR}/data/${VCF_HG38_PREFIX}_${chr_name}_variants.txt ${DIR}/data/postMIS_hg38_${chr_name}_variants.txt | wc -l | cut -d' ' -f1)
echo -e "${chr_name}\t${preMIS_only}\t${common}\t${postMIS_only}" > ${DIR}/data/preMIS_vs_postMIS_variant_counts_${chr_name}.txt

# Clean up
rm ${DIR}/data/${VCF_HG38_PREFIX}_${chr_name}.vcf.gz
