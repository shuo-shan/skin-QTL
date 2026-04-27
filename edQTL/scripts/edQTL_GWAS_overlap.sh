#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=2040]
#BSUB -q long
#BSUB -W 08:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

# set-up working directory
dir=/pi/manuel.garber-umw/human/skin/eQTLs/edQTL
cd ${dir}/output

# enlish edQTLs from pipeline output
cat ${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal | cut -d' ' -f1,6,9,11 | tr ' ' '\t' > edSite_edQTL_slope_padj.txt

# inner-join with GWAS SNP table
gwas_table=/pi/manuel.garber-umw/human/skin/eQTLs/GWAS_SNPs/

# fetch high LD tags (r2>0.8) of GWAS SNPs
cat /pi/manuel.garber-umw/human/skin/eQTLs/GWAS_SNPs/gwas_catalog_v1.0_all_associations_snps.txt | awk 'NR>1' | cut -f4 | sort | uniq > gwas_snps.txt
while read snp;do
  bash compile_high_LD_tags_of_snp.sh ${snp}
done < gwas_snps.txt

# this resulted in a file called ${dir}/output/gwas_snps_and_LD_tags.txt
cut -f1 gwas_snps_and_LD_tags.txt > temp1
cut -f2 gwas_snps_and_LD_tags.txt > temp2
cat temp1 temp2 | sort | uniq > gwas_snps_and_LD_tags_long.txt




# step 7. Filter out significant results to plot in R
cd ${dir}/output
cat ${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal | awk '$11<0.05' > foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal.sig
cat foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal.sig | cut -d' ' -f6 | sort | uniq > temp.snps
vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
module load bcftools/1.16
bcftools view -h ${vcfF} > temp.header.txt
bcftools view --include ID==@temp.snps ${vcfF} -Oz -o temp.snps.vcf.gz
bcftools view -H temp.snps.vcf.gz | cat temp.header.txt - > temp.snps.vcf.with.header
bcftools view temp.snps.vcf.with.header -Oz -o snps.vcf.gz
rm temp.header.txt temp.snps.vcf.with.header













