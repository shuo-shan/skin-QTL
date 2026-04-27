#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=300000]
#BSUB -q long
#BSUB -W 72:00
#BSUB -J plinkQC
#BSUB -e "./plinkQC_individual_%J.err"
#BSUB -o "./plinkQC_individual_%J.out"

### --------- in the cluster ----------- 
dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC
cd ${dir}
mkdir -p genodata
mkdir -p genodata/qc

### first compile all necessary documents
module load bcftools
module load plink/1.90b6.27

## Create PCA of merged dataset:
#cd ${dir}/genodata
#vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/merged_lowGPreset.vcf.gz
#plink --vcf ${vcf} --make-bed --out data_nosex
## assign sex to plink .fam file
#awk '{ $5 = 1; print }' data_nosex.fam > data_with_sex.fam
#mv data_with_sex.fam data_nosex.fam
#plink --bfile data_nosex --make-bed --out data
#rm test_nosex*
#plink --bfile data --pca --out data
#plink --bfile data --sexcheck --out data
#rm data_nosex*
#
## HapMap data and genotype PCA
#cd ${dir}/genodata/qc
#wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/HapMap_ID2Pop.txt -O HapMap_ID2Pop.txt
#wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/HapMap_PopColors.txt -O HapMap_PopColors.txt
#wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/data.HapMapIII.eigenvec -O data.HapMapIII.eigenvec


### --------- Rscript  ----------- 
module load r/4.2.2
Rscript run_plinkQC_per_individual_QC.R
