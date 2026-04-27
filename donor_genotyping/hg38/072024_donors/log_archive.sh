
dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors/selected_with_inDropData
cd ${dir}
module load bcftools

# modify GT tag to ./. if lowquality SNP for each donor's VCF file`
cd /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors/organized; module load bcftools; for f in *.vcf.gz;do bcftools +setGT ${f} -- -t q -i 'FILTER="LOWQUALSNP"' -n './.' | bcftools view -Oz -o modified_${f}; bcftools index modified_${f}; echo ${f}; done

cd merged
# retrieve autosomes
regions=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22
bcftools view -r ${regions} merged.vcf.gz -Oz -o merged_autosomes.vcf.gz
bcftools index merged_autosomes.vcf.gz

# retrieve biallelic SNPs
bcftools view -v snps -m2 -M2  merged_autosomes.vcf.gz -Oz -o merged_autosomes_biallelicSNPs.vcf.gz
bcftools index merged_autosomes_biallelicSNPs.vcf.gz

# query just reQTLs or GWAS SNPs of interest
bcftools view -i 'ID=@query_snps/all_candidate_reQTLs_052024_also_skinDis_melanoma_autoImmuneDis_snps.txt' merged_autosomes_biallelicSNPs.vcf.gz -Oz -o merged_autosomes_biallelicSNPs_all_candidate_reQTLs_052024_also_skinDis_melanoma_autoImmuneDis_snps.vcf.gz

bcftools view -i 'ID=@skinDis_melanoma_autoImmuneDis_snps.txt' merged_autosomes_biallelicSNPs.vcf.gz -Oz -o merged_autosomes_biallelicSNPs_skinDis_melanoma_autoImmuneDis_snps.vcf.gz

# get genotype
bcftools query -f '%CHROM %POS %ID %REF %ALT[ %GT]\n' merged_autosomes_biallelicSNPs.vcf.gz > merged_autosomes_biallelicSNPs.genotype.txt


# compare low cov and 10x cov missing rates for skin disease GWAS SNPs and reQTLs.
cp /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.vcf.gz ./round1_merged.vcf.gz

cd /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors/merged; module load bcftools; bcftools view -i 'ID=@skinDis_melanoma_autoImmuneDis_snps.txt' round1_merged.vcf.gz -Oz -o round1_merged_skinDis_melanoma_autoImmuneDis_snps.vcf.gz; bcftools index round1_merged_skinDis_melanoma_autoImmuneDis_snps.vcf.gz; bcftools query -f '%CHROM %POS %ID %REF %ALT[ %GT]\n' round1_merged_skinDis_melanoma_autoImmuneDis_snps.vcf.gz > round1_merged_skinDis_melanoma_autoImmuneDis_snps.txt
cd /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors/merged; module load bcftools; bcftools view -i 'ID=@skinDis_melanoma_autoImmuneDis_snps.txt' modified_merged.vcf.gz -Oz -o merged_skinDis_melanoma_autoImmuneDis_snps.vcf.gz; bcftools index merged_skinDis_melanoma_autoImmuneDis_snps.vcf.gz; bcftools query -f '%CHROM %POS %ID %REF %ALT[ %GT]\n' merged_skinDis_melanoma_autoImmuneDis_snps.vcf.gz > merged_skinDis_melanoma_autoImmuneDis_snps.txt
cd /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors/merged; module load bcftools; bcftools view -i 'ID=@all_candidate_reQTLs_052024.txt' modified_merged.vcf.gz -Oz -o merged_all_candidate_reQTLs_052024.vcf.gz; bcftools index merged_all_candidate_reQTLs_052024.vcf.gz; bcftools query -f '%CHROM %POS %ID %REF %ALT[ %GT]\n' merged_all_candidate_reQTLs_052024.vcf.gz > merged_all_candidate_reQTLs_052024.txt
cd /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors/merged; module load bcftools; bcftools view -i 'ID=@all_candidate_reQTLs_052024.txt' round1_merged.vcf.gz -Oz -o round1_merged_all_candidate_reQTLs_052024.vcf.gz; bcftools index round1_merged_all_candidate_reQTLs_052024.vcf.gz; bcftools query -f '%CHROM %POS %ID %REF %ALT[ %GT]\n'  round1_merged_all_candidate_reQTLs_052024.vcf.gz > round1_merged_all_candidate_reQTLs_052024.txt
