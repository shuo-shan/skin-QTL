#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=204000]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

cd /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/HWE
##vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
##module load vcftools
##vcftools --gzvcf ${vcfF} --hardy
##mv out.hwe /bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.hwe
#
#vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.vcf.gz
## only select European donors
#module load bcftools
#bcftools view -S european_donors.txt -o european_variants.vcf.gz -O z ${vcfF}
## only select for biallelic snps
#bcftools view -v snps -m2 -M2 -o european_biallelic_snps.vcf.gz -O z european_variants.vcf.gz
#module load vcftools
#vcftools --gzvcf european_biallelic_snps.vcf.gz --hardy
#mv out.hwe bd841628-fcc2-487a-8460-f5428237f0c9.merged.european_biallelic_snps.hwe
#cat bd841628-fcc2-487a-8460-f5428237f0c9.merged.european_biallelic_snps.hwe | awk 'NR==1{print $6}NR>1{ if ($5!="-nan") print $6}' > european_biallelic_snps_HWEpval.txt
#conda activate fastQTL
#Rscript make_histogram.R


## get variants in HWE and then perform pruning
vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/HWE/european_biallelic_snps.vcf.gz
module load bcftools
bcftools view ${vcfF} | bcftools +fill-tags -- -t HWE | bcftools filter -i 'INFO/HWE>0.000001' -O z -o variants_in_HWE.vcf.gz 
module load plink
plink --vcf variants_in_HWE.vcf.gz --make-bed --out data
plink --bfile data --indep-pairwise 100 5 0.2 --out pruned_data
plink --bfile data --extract pruned_data.prune.in --make-bed --out final_pruned_data
plink --bfile final_pruned_data --recode vcf --out final_pruned_data
bcftools view final_pruned_data.vcf -O z -o variants_in_HWE_pruned.vcf.gz
# clean up
rm data.bed data.bim data.fam data.log data.nosex final_pruned_data*

## separate variants by MAF tier
vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/HWE/variants_in_HWE_pruned.vcf.gz
module load bcftools
module load plink
plink --vcf ${vcfF} --make-bed --out data
plink --bfile data --freq --out freq_output
# 0-0.1 (--maf is inclusive, --max-maf is exclusive)
plink --bfile data --max-maf 0.1 --make-bed --out temp_plink_bed
plink --bfile temp_plink_bed --recode vcf --out variants_in_HWE_pruned_maftier
bcftools view variants_in_HWE_pruned_maftier.vcf -O z -o variants_in_HWE_pruned_maf0to0.1.vcf.gz
rm temp_plink_bed* variants_in_HWE_pruned_maftier*
# 0.1-0.2
plink --bfile data --maf 0.1 --max-maf 0.2 --make-bed --out temp_plink_bed
plink --bfile temp_plink_bed --recode vcf --out variants_in_HWE_pruned_maftier
bcftools view variants_in_HWE_pruned_maftier.vcf -O z -o variants_in_HWE_pruned_maf0.1to0.2.vcf.gz
rm temp_plink_bed* variants_in_HWE_pruned_maftier*
# 0.2-0.3
plink --bfile data --maf 0.2 --max-maf 0.3 --make-bed --out temp_plink_bed
plink --bfile temp_plink_bed --recode vcf --out variants_in_HWE_pruned_maftier
bcftools view variants_in_HWE_pruned_maftier.vcf -O z -o variants_in_HWE_pruned_maf0.2to0.3.vcf.gz
rm temp_plink_bed* variants_in_HWE_pruned_maftier*
# 0.3-0.4
plink --bfile data --maf 0.3 --max-maf 0.4 --make-bed --out temp_plink_bed
plink --bfile temp_plink_bed --recode vcf --out variants_in_HWE_pruned_maftier
bcftools view variants_in_HWE_pruned_maftier.vcf -O z -o variants_in_HWE_pruned_maf0.3to0.4.vcf.gz
rm temp_plink_bed* variants_in_HWE_pruned_maftier*
# 0.4-0.5
plink --bfile data --maf 0.4 --make-bed --out temp_plink_bed
plink --bfile temp_plink_bed --recode vcf --out variants_in_HWE_pruned_maftier
bcftools view variants_in_HWE_pruned_maftier.vcf -O z -o variants_in_HWE_pruned_maf0.4to0.5.vcf.gz
rm temp_plink_bed* variants_in_HWE_pruned_maftier*
# clean-up
rm data*

## overlap with
for f in variants_in_HWE_pruned_maf*;do
	vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/HWE/${f}
	prefix=$(echo ${f} | )
	res=modeling_results_featureSelectedModel_phenotypeRankNormCPM.txt
	bcftools view -H ${f} | awk '{print $3}' | sort | uniq > 

done




### archive
#fName=snps_used_in_DREG_KRT
#bcftools view --include ID==@${fName}.txt ${vcfF} | bcftools +fill-tags -- -t AN,AC_Hom,AC_Het,ExcHet,HWE | bcftools query -f '%ID\t%REF\t%ALT\t%AN\t%AC_Het\t%AC_Hom\t%HWE\t%ExcHet\n' > ${fName}_bcftools_query.txt 
#sleep 10
#fName=snps_used_in_DREG_MEL
#bcftools view --include ID==@${fName}.txt ${vcfF} | bcftools +fill-tags -- -t AN,AC_Hom,AC_Het,ExcHet,HWE | bcftools query -f '%ID\t%REF\t%ALT\t%AN\t%AC_Het\t%AC_Hom\t%HWE\t%ExcHet\n' > ${fName}_bcftools_query.txt
#sleep 10
#fName=snps_used_in_DREG_FRB
#bcftools view --include ID==@${fName}.txt ${vcfF} | bcftools +fill-tags -- -t AN,AC_Hom,AC_Het,ExcHet,HWE | bcftools query -f '%ID\t%REF\t%ALT\t%AN\t%AC_Het\t%AC_Hom\t%HWE\t%ExcHet\n' > ${fName}_bcftools_query.txt
## it is necessary to add another column of homozygous reference allele count (AN-
#
