module load bcftools/1.4.1
bcftools view --min-ac 1 renamed_impute-vcf-merged.vcf.bgz --threads 8 -O z >filt_renamed_impute-vcf-merged.vcf.gz
bcftools +fill-tags filt_renamed_impute-vcf-merged.vcf.gz -o filt_renamed_impute-vcf-merged_with_AF.vcf.gz -- -t AF
bgzip filt_renamed_impute-vcf-merged_with_AF.vcf
tabix -p vcf filt_renamed_impute-vcf-merged_with_AF.vcf.gz
