module load bcftools
vcf="/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz"
snp=rs8096411
bcftools view -i "ID==\"${snp}\"" ${vcf} | tail -2 > filtered_${snp}.txt
