#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=100000]
#BSUB -q long
#BSUB -W 24:00
#BSUB -e "./vcf2bed.%J%I.err"
#BSUB -o "./vcf2bed.%J%I.out"

dir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG
celltype=FRB
jobname=${celltype}
genotype_table=/pi/manuel.garber-umw/human/skin/eQTLs/fastQTL/GWAS_SNPs/SNPs_in_RE/snps_in_gene_proximal_${celltype}_in_regulatory_region.vcf.gz

module load bcftools/1.9
module load condas/2018-05-11
source activate fastQTL

# obtain the genotype and info from vcf file 
bcftools query -f "%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT{0}[\t%GT]\n" -H ${genotype_table} | head -1 > header
cat header | tr '\t' '\n' | cut -d']' -f2 | tr '\n' '\t'  | sed '$s/\t$/\n/' > header2
cat header2 | sed 's/POS\tPOS/START\tEND/g'  > header3
bcftools query -f "%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT{0}[\t%GT]\n" ${genotype_table} -o temp.filtered.genotype.vcf
cat temp.filtered.genotype.vcf | awk '{OFS="\t"}{print $1,$2,$2+1,$4,$5,$6}' > temp1
cat temp.filtered.genotype.vcf | cut -f7- > temp2
paste temp1 temp2 > temp3
cat temp3 > temp.filtered.genotype.bed
cat header3 > temp.header
rm header header2 header3 temp.filtered.genotype.vcf temp1 temp2 temp3
echo "there are " $(wc -l temp.filtered.genotype.bed) " SNPs in functional emement regions".
# 3934250 SNPs in functional emement regions
# 3910974 SNPs in KRT functional emement regions
# 3936539 SNPs in FRB functional emement regions 

### further filtering (3934250 snps -->  524406 snps) 
# this ensures that all the SNPs we are looking at have at least 3 donors in each of the 3 genotypes.
# split bed file into 50K entries per file.
split -l 50000 temp.filtered.genotype.bed temp.split.
sleep 10
rm commands.txt
for f in temp.split.*;do
  echo "date;Rscript /pi/manuel.garber-umw/human/skin/eQTLs/fastQTL/fastQTL/GWAS_SNP_with_LD/further_filtering_SNPs_ver2.R ${dir} ${f} temp.header; date" >> commands.txt
done
while read c; do
  echo ${c} | bsub -W 8:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long -e "./log/further_filtering.%J%I.err" -o "./log/further_filtering.%J%I.out"
done < commands.txt
while [[ $(bjobs | grep ${jobname} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${jobname} | wc -l) "jobs remaining"; sleep 5; done
# all outputs are named new.xxxx. without a header. join them and then add a header.
cat temp.header new.temp.split.* > filtered.genotype.bed
mv filtered.genotype.bed ${celltype}_filtered.genotype.bed
rm temp.header temp.split.* new.temp.split.* temp.filtered.genotype.bed commands.txt
echo "there are " $(wc -l ${celltype}_filtered.genotype.bed) " SNPs in functional emement regions".
# this results in a filtered.genotype.bed that contain SNPs with all 3 genotypes, and at least 3 donors per genotype.
# 521,329 filtered KRT SNPs
# 524,406 filtered MEL SNPs
# 524,616 filtered FRB SNPs

