#!/bin/bash
# written by Crystal Shan 06/2022, modified on 10/2022
# goal: for a given SNP, compile a data matrix, build linear model, then output modeling results.
# output: a folder containing this snp and relevant plots and genes.

### inputs
snp=$1
covariates=$2 #a1b4c7 a3b3c9
dir=$3

### 1. set-up
celltype=FRB
genotype_table=/pi/manuel.garber-umw/human/skin/eQTLs/edQTL/output/genotype.bed
datadir=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis
PBS_table=${datadir}/CPM_expressedGenes_PBS.txt
IFN_table=${datadir}/CPM_expressedGenes_IFN.txt
expressed_genes=${datadir}/expressedGenes.txt
conda activate fastQTL
### 2. create and go to temp folder
cd ${dir}/modelingResult_${covariates}
mkdir ${dir}/modelingResult_${covariates}/temp_${snp}
cd ${dir}/modelingResult_${covariates}/temp_${snp}

### 3. create genotype file of the SNP
cat ${genotype_table} | head -1 |  tr '\t' '\n' | sed 's/\:GT//g' > temp.header
cat ${genotype_table} | grep -w ${snp}  | tr '\t' '\n' | paste temp.header - > temp.${snp}
chr=$(grep -w "CHROM" temp.${snp} | cut -f2)
start_pos=$(grep -w "START" temp.${snp} | cut -f2)
echo -e "donor\tgenotype" > genotype.txt
cat temp.${snp} | awk 'NR>6' | sed 's/0\/0/0/g' | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' | sed 's/\.\/\./NA/g' | grep -v NA >> genotype.txt
rm temp.header temp.${snp}
echo "check point 1: made snp genotype file"; date;

### 4. fetch nearby expressed genes of the SNP
cat ${expressed_genes} | grep -w ${chr} | awk -v p=${start_pos} '{if ($2>=p-500000 && $2<p+500000) print $5}' > expressed_nearby_genes.txt
while read g;do
	echo -e "donor\tPBS" > PBS_${g}.txt
	echo -e "donor\tIFN" > IFN_${g}.txt
	head -1 ${PBS_table} | cut -f2- | tr '\t' '\n' > PBS_${g}_header
	grep -w ${g} ${PBS_table} | cut -f2- | tr '\t' '\n' | paste PBS_${g}_header - | sort -k1 >> PBS_${g}.txt
	head -1 ${IFN_table} | cut -f2- | tr '\t' '\n' > IFN_${g}_header
	grep -w ${g} ${IFN_table} | cut -f2- | tr '\t' '\n' | paste IFN_${g}_header - | sort -k1 >> IFN_${g}.txt
	rm PBS_${g}_header IFN_${g}_header
done < expressed_nearby_genes.txt
echo "check point 2: fetched nearby genes rank Norm scores"; date;

## 5. build linear regression model for snp-gene pair
while read g;do
	echo "working on SNP:gene pair "${snp}":"${g}; date;
        this_dir=${dir}/modelingResult_${covariates}/temp_${snp}
	PBS_CPM=${this_dir}/PBS_${g}.txt
	IFN_CPM=${this_dir}/IFN_${g}.txt
	genotype=${this_dir}/genotype.txt
	covariates_matrix=${dir}/covariates_${covariates}.txt
	Rscript ${dir}/build_linear_model_log2FC.R ${dir} ${snp} ${g} ${covariates} ${PBS_CPM} ${IFN_CPM} ${genotype} ${covariates_matrix} reQTL
done < ${dir}/modelingResult_${covariates}/temp_${snp}/expressed_nearby_genes.txt
cd ${dir}
rm -r ${dir}/modelingResult_${covariates}/temp_${snp}
echo "check point 4: built models for each snp-gene pair of snp:"${snp}; date;

