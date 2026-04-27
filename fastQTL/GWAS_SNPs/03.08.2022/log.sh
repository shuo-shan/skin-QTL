# script written by Crystal SHan 03/2022
# goal: filter donor vcf file with GWAS SNPs

# In this folder:
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/03.08.2022
cat skin_disease_snps.txt | sed 's/chr//g' | sort | uniq > skin_disease_GWAS_SNPs.txt
### get all GWAS SNPs for skin disease
module load condas/2018-05-11
source activate fastQTL
module load bcftools/1.9  
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
# get header
bcftools view -h ${genotypeF} > header
# get body
bcftools view -H ${genotypeF} > body
# query body
dir=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/03.08.2022
gwas_snps=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/03.08.2022/skin_disease_GWAS_SNPs.txt
bash /nl/umw_manuel_garber/sshan/scripts/function_rapid_fgrep.sh ${gwas_snps} ${dir}/body temp.filtered_skin_disease_GWAS_SNPs.vcf name1 100000
cat header temp.filtered_skin_disease_GWAS_SNPs.vcf > filtered_skin_disease_GWAS_SNPs.vcf
rm temp.filtered_skin_disease_GWAS_SNPs.vcf header body; rm -r rapid_fgrep_temp/
### turn VCF file into bed file
bcftools query -f "%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT{0}[\t%GT]\n" -H filtered_skin_disease_GWAS_SNPs.vcf | head -1 > header
cat header | tr '\t' '\n' | cut -d']' -f2 | tr '\n' '\t'  | sed '$s/\t$/\n/' > header2
cat header2 | sed 's/POS\tPOS/START\tEND/g'  > header3
bcftools query -f "%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT{0}[\t%GT]\n" filtered_skin_disease_GWAS_SNPs.vcf -o temp.filtered.genotype.vcf
cat temp.filtered.genotype.vcf | awk '{OFS="\t"}{print $1,$2,$2+1,$4,$5,$6}' > temp1
cat temp.filtered.genotype.vcf | cut -f7- > temp2
paste temp1 temp2 > temp3
cat header3 temp3 > temp.filtered.genotype.bed
rm header header2 header3 temp.filtered.genotype.vcf temp1 temp2 temp3
## further filtering ( snps --> 29,951 snps) 
#Rscript /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/GWAS_SNP_with_LD/further_filtering_SNPs_ver2.R ${dir} temp.filtered.genotype.bed
#conda deactivate
## ^ results in a filtered.genotype.bed: SNPs with at least 3 genotype groups (HomoRef, Het, HomoAlt) and at least 3 donors in each group
#rm temp.filtered.genotype.bed

############################# 10/21/2021
### Get SNPs in high LD
### Get all variants in high LD (r2 > 0.6) for all GWAS SNPs within a range of 1000kb
# script: /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/plink/log.sh
cd /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/plink
pref=snps_plus_high_ld
gwas_snps=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/skin_disease_GWAS_SNPs.txt
plink --file 37donors --r2 --ld-snp-list ${gwas_snps} --ld-window-kb 1000 --ld-window-r2 0.6 --out ${pref}
# SNPs with high LD info are stored in: 37donors_snps_plus_high_ld.ld
### organize GWAS SNPs and all those in high LD into a single column
cat snps_plus_high_ld.ld |  awk '{OFS=";"}NR>1{print $3,$6}' | tr ';' '\n' | sort | uniq > 37donors_snps_plus_high_ld_lst.txt
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/03.08.2022
ln -s /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/plink/37donors_snps_plus_high_ld_lst.txt
###
temp=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/plink/37donors_snps_plus_high_ld_lst.txt
cat $temp | awk '{print $0"\t"}' > temp
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
bcftools view -h ${genotypeF} > header
bcftools view -H ${genotypeF} > body
split -l 100000 body temp.split.
rm commands.txt
for f in temp.split.*;do
  echo "grep -w -f ${temp} ${f} > grepped.${f}" >> commands.txt
done
while read c;do
  echo ${c} | bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long -e "./log/%J%I.err" -o "./log/%J%I.out"
done < commands.txt
cat grepped.temp.split.* > tempbody
cat header tempbody > 37donors_skin_disease_GWAS_SNPs_plus_high_ld.vcf
bgzip 37donors_skin_disease_GWAS_SNPs_plus_high_ld.vcf && tabix -p vcf 37donors_skin_disease_GWAS_SNPs_plus_high_ld.vcf.gz 

rm grepped.temp.split.*
rm temp.split.*
rm temp header body tempbody filtered_skin_disease_GWAS_SNPs.vcf
rm commands.txt

############################# 01/24/2022
### get Minor Allele Frequency from VCF file
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs
module load condas/2018-05-11
source activate sshan_isoform
module load bcftools/1.9
f=filtered_skin_disease_GWAS_SNPs.vcf
bcftools view -H ${f} | head
bcftools view -H ${f} | cut -f8 | sed 's/.*AF=//g' > temp.maf
bcftools view -H ${f} | cut -f3,1,2  > temp.snp
paste temp.snp temp.maf > filtered_skin_disease_GWAS_SNPs_MAF.txt
rm temp.maf temp.snp

### get Minor Allele Frequency from VCF file filtered from master VCF
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


