#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 250:00
#BSUB -e "./pipeline.%J%I.err"
#BSUB -o "./pipeline.%J%I.out"

module load bedtools/2.30.0
celltype=KRT
snp=rs6961406
gene=SNHG26
jobname=${celltype}
dir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/gene_SNP_of_interest/${gene}_${snp}
mkdir ${dir}
cd ${dir}

####### 1. compile a bed file for the SNP genotype file.
echo ${snp} > snp.txt
master_vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
module load bcftools/1.16
source activate fastQTL
bcftools view --include ID==@snp.txt ${master_vcf} -Oz -o ${dir}/snp.vcf.gz
genotype_table=${dir}/snp.vcf.gz
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
cat temp.header > snp.bed
cat temp.filtered.genotype.bed >> snp.bed
mv temp.filtered.genotype.bed temp.header

####### 2. Compile log2FC covariate matrix based on the desired number of genotypePCs and latent phenotype variables.
for n_PCs in 10;do
    for n_latentVar in 10;do
      mkdir ${dir}/modelingResult_a${n_PCs}b${n_latentVar}
      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} CPM_PBS
      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} CPM_IFN
      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} Log2FCwithDummy10
      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} rankNormCPM_PBS
      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} rankNormCPM_IFN
      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} rankNormLog2FCwithDummy10
    done
done

######### 3. For every SNP of interest and gene of interest, construct the model
### 3.1 set-up
n_PCs=10
n_latentVar=10
covariates=a10b10
genotype_table=${dir}/snp.bed
datadir=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis
PBS_table=${datadir}/CPM_expressedGenes_PBS.txt
IFN_table=${datadir}/CPM_expressedGenes_IFN.txt
source activate fastQTL

### 3.2 create and go to temp folder
cd ${dir}/modelingResult_${covariates}
mkdir ${dir}/modelingResult_${covariates}/temp_${snp}
cd ${dir}/modelingResult_${covariates}/temp_${snp}

### 3.3 create genotype file of the SNP
cat ${genotype_table} | head -1 |  tr '\t' '\n' | sed 's/\:GT//g' > temp.header
cat ${genotype_table} | grep -w ${snp}  | tr '\t' '\n' | paste temp.header - > temp.${snp}
chr=$(grep -w "CHROM" temp.${snp} | cut -f2)
start_pos=$(grep -w "START" temp.${snp} | cut -f2)
echo -e "donor\tgenotype" > genotype.txt
cat temp.${snp} | awk 'NR>6' | sed 's/0\/0/0/g' | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' | sed 's/\.\/\./NA/g' | grep -v NA >> genotype.txt
rm temp.header temp.${snp}
echo "check point 1: made snp genotype file"; date;

### 3.4 fetch gene expression for the gene of interest
echo ${gene} > gene.txt
while read g;do
        echo -e "donor\tPBS" > PBS_${g}.txt
        head -1 ${PBS_table} | tr '\t' '\n' | awk 'NF' > PBS_${g}_header
        awk -v gene=${g} '$1==gene' ${PBS_table} | cut -f2- | tr '\t' '\n' | awk 'NF' | paste PBS_${g}_header - | sort -k1 >> PBS_${g}.txt

        echo -e "donor\tIFN" > IFN_${g}.txt
        head -1 ${IFN_table} | tr '\t' '\n' | awk 'NF' > IFN_${g}_header
        awk -v gene=${g} '$1==gene' ${IFN_table} | cut -f2- | tr '\t' '\n' | awk 'NF' | paste IFN_${g}_header - | sort -k1 >> IFN_${g}.txt
        rm PBS_${g}_header IFN_${g}_header
done < gene.txt
echo "check point 2: fetched gene of interest CPM values in PBS and IFN"; date;

### 3.5 build linear regression model for snp-gene pair
echo "working on SNP:gene pair "${snp}":"${g}; date;
cd ${dir}/modelingResult_${covariates}/temp_${snp}
this_dir=${dir}/modelingResult_${covariates}/temp_${snp}
PBS_CPM=${this_dir}/PBS_${gene}.txt
IFN_CPM=${this_dir}/IFN_${gene}.txt
genotype=${this_dir}/genotype.txt
Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/build_linear_model_with_rankNorm_and_permutation_09072023.R ${dir} ${snp} ${gene} ${covariates} ${PBS_CPM} ${IFN_CPM} ${genotype}
cd ${dir}/modelingResult_${covariates}

### 3.6 clean up temporary files/folders
rm ${snp}*.txt
rm -r ${dir}/modelingResult_${covariates}/temp_${snp}
cd ${dir}
echo "check point 4: built models for each snp-gene pair of snp:"${snp}; date;


####### 4. Make the plots
# go to /pi/manuel.garber-umw/human/skin/eQTLs/DREG/SNP_of_interest/make_DREG_plots_for_website.R





















