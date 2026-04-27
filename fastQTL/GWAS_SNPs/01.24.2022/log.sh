#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=18000]
#BSUB -q long
#BSUB -W 21:00
#BSUB -o “./%J%I.out” 
#BSUB -e “./%J%I.err”
### comment

cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/01.24.2022
module load condas/2018-05-11
source activate sshan_isoform
module load bcftools/1.9
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.vcf.gz
# get header
bcftools view -h ${genotypeF} > header
# get body
bcftools view -H ${genotypeF} > body
# query body
split -d -l 100000 body temp.split.
for f in temp.split.*;do
  echo "grep -f skin_disease_GWAS_SNPs.txt ${f} > grepped.${f}" >> commands.txt
done
while read c;do
  echo ${c} |     bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long
done < commands.txt
cat grepped.temp.split.* > tempbody
cat header tempbody > filtered_skin_disease_GWAS_SNPs.vcf
bcftools +fill-tags filtered_skin_disease_GWAS_SNPs.vcf -- -t 'MAF' > filtered_skin_disease_GWAS_SNPs.tagged.vcf
rm grepped.* temp.split.* header body tempbody
#mv skin_disease_GWAS_SNPs.txt temp; cat temp | sed 's/\t//g' > skin_disease_GWAS_SNPs.txt # remove the tab at the end
#mv skin_disease_GWAS_SNPs.txt temp; cat temp | sed 's/ rs/rs/g' | sort | uniq > skin_disease_GWAS_SNPs.txt
f=filtered_skin_disease_GWAS_SNPs.tagged.vcf
bcftools view -H ${f} | head
bcftools view -H ${f} | cut -f8 | sed 's/.*AF=//g' > temp.maf
bcftools view -H ${f} | cut -f3,1,2  > temp.snp
paste temp.snp temp.maf > filtered_skin_disease_GWAS_SNPs_MAF.txt
rm temp.maf temp.snp
