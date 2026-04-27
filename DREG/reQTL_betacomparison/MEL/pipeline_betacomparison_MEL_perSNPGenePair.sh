#!/bin/bash
# written by Crystal Shan 06/2022, modified on 10/2022
# goal: for a given SNP, compile a data matrix, build linear model, then output modeling results.
# output: a folder containing this snp and relevant plots and genes.

### inputs
snp=$1
g=$2 # gene
dir=$3
outputName=$4

### 1. set-up
celltype=MEL
genotype_table=${dir}/snps_near_expressed_genes.bed
datadir=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis
PBS_table=${datadir}/CPM_expressedGenes_PBS.txt
IFN_table=${datadir}/CPM_expressedGenes_IFN.txt
expressed_genes=${datadir}/expressedGenes.txt
source activate fastQTL
### 2. create and go to temp folder
mkdir ${dir}/modelingResult/temp_${snp}_${g}
cd ${dir}/modelingResult/temp_${snp}_${g}

### 3. create genotype file of the SNP
cat ${genotype_table} | head -1 |  tr '\t' '\n' | sed 's/\:GT//g' > temp.header
cat ${genotype_table} | grep -w ${snp}  | tr '\t' '\n' | paste temp.header - > temp.${snp}
chr=$(grep -w "CHROM" temp.${snp} | cut -f2)
start_pos=$(grep -w "START" temp.${snp} | cut -f2)
echo -e "donor\tgenotype" > genotype.txt
cat temp.${snp} | awk 'NR>6' | sed 's/0\/0/0/g' | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' | sed 's/\.\/\./NA/g' | grep -v NA >> genotype.txt
rm temp.header temp.${snp}
echo "check point 1: made snp genotype file"; date;

### 4. fetch gene expression table
echo -e "donor\tPBS" > PBS_${g}.txt
head -1 ${PBS_table} | tr '\t' '\n' | awk 'NF' > PBS_${g}_header
awk -v gene=${g} '$1==gene' ${PBS_table} | cut -f2- | tr '\t' '\n' | awk 'NF' | paste PBS_${g}_header - | sort -k1 >> PBS_${g}.txt

echo -e "donor\tIFN" > IFN_${g}.txt
head -1 ${IFN_table} | tr '\t' '\n' | awk 'NF' > IFN_${g}_header
awk -v gene=${g} '$1==gene' ${IFN_table} | cut -f2- | tr '\t' '\n' | awk 'NF' | paste IFN_${g}_header - | sort -k1 >> IFN_${g}.txt
rm PBS_${g}_header IFN_${g}_header
echo "check point 2: fetched genes CPM values in PBS and IFN"; date;

## 5. build beta comparison model for snp-gene pair
echo "working on SNP:gene pair "${snp}":"${g}; date;
this_dir=${dir}/modelingResult/temp_${snp}_${g}
PBS_CPM=${this_dir}/PBS_${g}.txt
IFN_CPM=${this_dir}/IFN_${g}.txt
genotype=${this_dir}/genotype.txt
Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/build_betacomparison_model_with_permutation.R ${dir} ${snp} ${g} ${PBS_CPM} ${IFN_CPM} ${genotype}
cd ${dir}/modelingResult

# acquire lock before writing
output_file="${dir}/${outputName}.txt"
lock_file="${dir}/${outputName}.lock"
exec 9>"$lock_file" # this line opens the lock file for writing and associates it with file descriptor 9. lock file is created if it does not already exist.
flock 9 # apply an advisory lock on this file

# write to the output file
cat ${snp}"_"${g}.txt >> ${output_file}

# release the lock
flock -u 9

# clean up temporary files/folders
rm ${snp}"_"${g}.txt
rm -r ${dir}/modelingResult/temp_${snp}_${g}
cd ${dir}
echo "check point 4: built models for each snp-gene pair of:"${snp} ${g}; date;

