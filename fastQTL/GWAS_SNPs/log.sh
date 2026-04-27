# In Dropbox:
scp skin_disease_GWAS_SNPs_cleaned_list.txt 'ss65w@ghpcc06.umassrc.org://nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs'

# In this folder:
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs
cat skin_disease_GWAS_SNPs_cleaned_list.txt | cut -f1 | awk 'NR>1' | sort | uniq | sed 's/chr//g' | awk '{print $0"\t"}' > skin_disease_GWAS_SNPs.txt
# manually remove empty line
### get all GWAS SNPs for skin disease
module load condas/2018-05-11
source activate sshan_isoform
module load bcftools/1.9  
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.vcf.gz
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
bcftools view filtered_skin_disease_GWAS_SNPs.vcf -Oz -o filtered_skin_disease_GWAS_SNPs.vcf.gz
rm grepped.* temp.split.* header body tempbody
mv skin_disease_GWAS_SNPs.txt temp; cat temp | sed 's/\t//g' > skin_disease_GWAS_SNPs.txt # remove the tab at the end
mv skin_disease_GWAS_SNPs.txt temp; cat temp | sed 's/ rs/rs/g' | sort | uniq > skin_disease_GWAS_SNPs.txt

############################# 10/21/2021
### Get SNPs in high LD
### Get all variants in high LD (r2 > 0.6) for all GWAS SNPs within a range of 1000kb
# script: /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/plink/log.sh
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/GWAS_SNP_with_LD
temp=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/plink/28donors_snps_with_ld_lst.txt
cat $temp | awk '{print $0"\t"}' > temp
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.vcf.gz
lst=temp
bcftools view -h ${genotypeF} > header
bcftools view -H ${genotypeF} > body
split -d -l 100000 body temp.split.
for f in temp.split.*;do
  echo "grep -f ${lst} ${f} > grepped.${f}" >> commands.txt
done
while read c;do
  echo ${c} |     bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long
done < commands.txt
cat grepped.temp.split.* > tempbody
cat header tempbody > filtered_skin_disease_GWAS_SNPs.vcf
bgzip filtered_skin_disease_GWAS_SNPs.vcf && tabix -p vcf filtered_skin_disease_GWAS_SNPs.vcf.gz 

rm grepped.temp.split.*
rm temp.split.*
rm temp header body tempbody filtered_skin_disease_GWAS_SNPs.vcf
rm commands.txt

#remove sites with more than 0.8 linkeage disequilibrium
#module load bcftools/1.9
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged
#f=bd841628-fcc2-487a-8460-f5428237f0c9
#inf=${dir}/${f}.merged.filtered.2.16donors.vcf.gz
#bcftools +prune -l 0.8 -w 50 ${inf} -Oz -o ${dir}/${f}.16donors.filtered2.pruned.vcf.gz
#
## get Minor Allele Frequency from VCF file filtered from master VCF
### get all GWAS SNPs for skin disease
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
#bcftools view filtered_skin_disease_GWAS_SNPs.vcf -Oz -o filtered_skin_disease_GWAS_SNPs.vcf.gz
rm grepped.* temp.split.* header body tempbody
#mv skin_disease_GWAS_SNPs.txt temp; cat temp | sed 's/\t//g' > skin_disease_GWAS_SNPs.txt # remove the tab at the end
#mv skin_disease_GWAS_SNPs.txt temp; cat temp | sed 's/ rs/rs/g' | sort | uniq > skin_disease_GWAS_SNPs.txt
f=filtered_skin_disease_GWAS_SNPs.vcf
bcftools view -H ${f} | head
bcftools view -H ${f} | cut -f8 | sed 's/.*AF=//g' > temp.maf
bcftools view -H ${f} | cut -f3,1,2  > temp.snp
paste temp.snp temp.maf > filtered_skin_disease_GWAS_SNPs_MAF.txt
rm temp.maf temp.snp


############################# 05/10/2022
### prune out those in high LD to open region SNPs (r<0.8)
### discard sites with more than 0.8 linkeage disequilibrium
module load bcftools/1.9
dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged
f=bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged
inf=${dir}/${f}.vcf.gz
bcftools +prune -l 0.8 -w 50 ${inf} -Oz -o ${dir}/${f}.pruned.vcf.gz
bcftools index ${dir}/${f}.pruned.vcf.gz



 
