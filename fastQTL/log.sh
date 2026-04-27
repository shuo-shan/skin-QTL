#!/bin/bash
#BSUB -n 26
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=18000]
#BSUB -q long
#BSUB -W 121:00
### script to interface Gencove to retrieve, upload, and process genotyping data

### crystal shan 09/2021
### bsub -Is -q interactive -W 8:00 -n1 -R rusage[mem=450000] -R "span[hosts=1]" /bin/bash

############################
#### packages & env
module load condas/2018-05-11
source activate sshan_isoform
#pip install bx_python # to run collapse_annotation.py below.
#conda install -c anaconda scipy $ to run eqtl_prepare_expression.py
#pip3 install qtl
module load bcftools/1.9
module load bedops/2.4.14-x86_64
#docker pull broadinstitute/gtex_eqtl:V8 # pull docker image of GTEx eQTL pipeline
dir=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL
cd ${dir}
############################
#### subset vcf to 29 samples
module load tabix/0.2.6
f=bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
#head -1 counts.MEL_PBS.txt | tr '\t' '\n' | cut -d'_' -f1 | grep -v 'F47' > samples.txt # remove genotype failed sample F47
bcftools view --threads 32 --min-ac 3 -S samples.txt ${f} -Oz -o filtered.vcf.gz
tabix -p vcf filtered.vcf.gz
############################
#### output dosage-only vcf file
#bcftools query -f "%CHROM\t%POS\t%ID[\t%DS]\n" -H filtered.vcf.gz | head -1 > header
#cat header | tr '\t' '\n' | cut -d']' -f2 | sed 's/:DS//g' | tr '\n' '\t' > header2	
bcftools query -f "%CHROM\t%POS\t%ID[\t%DS]\n" filtered.vcf.gz -o filtered.dosage.vcf
cat header2 filtered.dosage.vcf > filtered.dosage.txt
rm filtered.dosage.vcf
#bcftools +dosage filtered.vcf.gz > filtered.imputed-dosage.txt
############################
#### output genotype-only vcf file
#bcftools query -f "%CHROM\t%POS\t%ID[\t%GT]\n" filtered.vcf.gz -o filtered.genotype.vcf
#cat header2 filtered.genotype.vcf > filtered.genotype.txt
#rm header2 filtered.genotype.vcf
############################
#### sample_partipant_lookup
head -1 counts.MEL_PBS.txt | tr '\t' '\n' | grep -v 'F47_PBS' > temp
paste temp samples.txt > sample_participant_lookup.txt
############################
#### vcf_chr_list
tabix -l filtered.vcf.gz > vcf_chr_list.txt
############################
#### download gene model & collapse
# python script obtained from GTEx: https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py
cd fastQTL
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz
gunzip gencode.v34.annotation.gtf.gz
#python3 collapse_annotation.py gencode.v34.annotation.gtf gencode.v34.GRCh38.genes.gtf
############################
#### convert gene model to bed file
module load bedops/2.4.14-x86_64
gtf2bed < gencode.v34.annotation.gtf > gencode.v34.annotation.gtf.bed
cat gencode.v34.annotation.gtf.bed | awk '{if ($8=="gene") print $0}' | cut -f10 | sed 's/.*gene_name //' | sed 's/; .*//' | sed -e 's/^"//' -e 's/"$//' > part1
cat gencode.v34.annotation.gtf.bed | awk '{if ($8=="gene") print $1"\t"$2"\t"$3"\t"}' > part2
paste part1 part2 > gencode.v34.annotation.bed
rm part1 part2 gencode.v34.annotation.gtf.bed
############################
#### generate phenotype matrices
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
anno=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/gencode.v34.annotation.bed
# 1. log2(cpm+1) data
ln -s /nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/tmm_normalized_log2cpm.MEL_PBS.bed
bgzip tmm_normalized_log2cpm.MEL_PBS.bed && tabix -p bed tmm_normalized_log2cpm.MEL_PBS.bed.gz
# 2. baseline cpm data
ln -s /nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/tmm_normalized_cpm.MEL_PBS.bed
bgzip tmm_normalized_cpm.MEL_PBS.bed && tabix -p bed tmm_normalized_cpm.MEL_PBS.bed.gz
# 3. log2FC data
ln -s /nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/log2FC_p0.05.MEL.bed
bgzip log2FC_p0.05.MEL.bed && tabix -p bed log2FC_p0.05.MEL.bed.gz
# 4. IFNg stimulated cpm data
ln -s /nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/tmm_normalized_cpm.MEL_IFN.bed
bgzip tmm_normalized_cpm.MEL_IFN.bed && tabix -p bed tmm_normalized_cpm.MEL_IFN.bed.gz
# 5. FC data
ln -s /nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/FC_p0.05.MEL.bed
bgzip FC_p0.05.MEL.bed && tabix -p bed FC_p0.05.MEL.bed.gz

############################
#### install & run PEER to calculate covariate
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL
module load condas/2018-05-11
conda install -c bioconda r-peer
conda create --name peer python=3.5.4
conda activate peer
conda config --add channels conda-forge
conda config --add channels r
conda config --add channels bioconda
conda install zlib=1.2.8
conda install r=3.4.1
conda install r-peer

#### run PEER
module load condas/2018-05-11
source activate peer
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
Rscript run_PEER.R tmm_normalized_log2cpm.MEL_PBS.bed.gz tmm_normalized_log2cpm.MEL_PBS 15
Rscript run_PEER.R tmm_normalized_cpm.MEL_PBS.bed.gz tmm_normalized_cpm.MEL_PBS 15
Rscript run_PEER.R log2FC_p0.05.MEL.bed.gz log2FC_p0.05.MEL 15
Rscript run_PEER.R tmm_normalized_cpm.MEL_IFN.bed.gz tmm_normalized_cpm.MEL_IFN 15
Rscript run_PEER.R FC_p0.05.MEL.bed.gz FC_p0.05.MEL 15

############################
#### combine genotype top 5 PCs and peer factors into covaraite file
#### genotype PCA was done in /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/genotype_PCA.R
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
cat genotype_dosage_topPCs.txt log2FC_p0.05.MEL.PEER_covariates.txt > log2FC_p0.05_MEL_covariates.txt # manually tweak the column name
cat genotype_dosage_topPCs.txt tmm_normalized_cpm.MEL_PBS.PEER_covariates.txt > tmm_normalized_cpm_MEL_PBS_covariates.txt
cat genotype_dosage_topPCs.txt tmm_normalized_log2cpm.MEL_PBS.PEER_covariates.txt > tmm_normalized_log2cpm_MEL_PBS_covariates.txt
cat genotype_dosage_topPCs.txt tmm_normalized_cpm.MEL_IFN.PEER_covariates.txt > tmm_normalized_cpm.MEL_IFN_covariates.txt
cat genotype_dosage_topPCs.txt FC_p0.05.MEL.PEER_covariates.txt > FC_p0.05.MEL_covariates.txt

##############################
#### install FastQTL
module load condas/2018-05-11
source activate sshan_isoform
module load boost/1.75.0
module load gsl/2.6
module load g++/8.1.0
module load binutils/2.37
## (B.1) INSTALL MATH R LIBRARY
#(1) Download R source code: 		wget http://cran.r-project.org/src/base/R-3/R-3.2.0.tar.gz
#(2) Unzip R source code: 			tar xzvf R-3.2.0.tar.gz
#(3) Go to R source code folder: 	cd R-3.2.0
#(4) Configure Makefile: 			./configure
#(5) Go to R math library folder: 	cd src/nmath/standalone
#(6) Compile the code: 				make
#(7) Go 2 folder backward:			cd ../..
#(8) Save the current path:			RMATH=$(pwd)
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL
make cleanall && make RMATH=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/R-3.2.0/src # created a binary file in FastQTL/bin folder
export PATH=$PATH:/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin # add this bin to path

################################
#### run FastQTL 
export PATH=${PATH}:/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.vcf.gz
phenotypeF1=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_cpm.MEL_PBS.bed.gz
phenotypeF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm.MEL_PBS.bed.gz
phenotypeF3=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log2FC_p0.05.MEL.bed.gz
phenotypeF4=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_cpm.MEL_IFN.bed.gz
phenotypeF5=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FC_p0.05.MEL.bed.gz

covF1=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_cpm_MEL_PBS_covariates.txt
covF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_MEL_PBS_covariates.txt
covF3=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log2FC_p0.05_MEL_covariates.txt
covF4=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_cpm.MEL_IFN_covariates.txt
covF5=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FC_p0.05.MEL_covariates.txt
# default simple
fastQTL --vcf ${genotypeF} --bed ${phenotypeF1} --region chr1:80000-1000000 --out nominals.default.test.txt.gz
# change phenotype data to be normally distributed in N(0,1)
fastQTL --vcf ${genotypeF} --bed ${phenotypeF1} --region chr1:80000-1000000 --out nominals.quantile.txt.gz --normal
# add covariate
fastQTL --vcf ${genotypeF} --bed ${phenotypeF1} --region chr1:80000-1000000 --out nominals.quantile.txt.gz --normal --cov ${covF1}
# permutation-based testing
fastQTL --vcf ${genotypeF} --bed ${phenotypeF1} --region chr1:80000-1000000 --permute 1000 --out permutations.quantile.txt.gz --normal --cov ${covF1}

###########################
# TMM-normalized log2 CPM from Melanocyte-PBS samples
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF2} --permute 1000 --normal --cov ${covF2} --out permutations_tmm_normalized_log2cpm.MEL_PBS.txt.gz --commands 40 commands.40.txt
while read c; do
     echo "${c}" |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10000] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.out" |\
     echo "submitted job"
done < commands.40.txt
# make sure to rename files

###########################
# log2 fold-change to detect responseQTL
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF3} --permute 1000 --normal --cov ${covF3} --out permutations_log2FC_p0.05.MEL.txt.gz --commands 40 commands.log2FC.40.txt
while read c; do
     echo "${c}" |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.log2FC.out" |\
     echo "submitted job"
done < commands.log2FC.40.txt

###########################
# TMM-normalized CPM to detect eQTL
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF1} --permute 1000 --normal --cov ${covF1} --out permutations_tmm_normalized_cpm.MEL_PBS.txt.gz --commands 40 commands.cpm.40.txt
while read c; do
     echo ${c} |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.cpm.out" |\
     echo "submitted job"
done < commands.cpm.40.txt

###########################
# log2 cpm, without covariates
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF1} --permute 1000 --normal --out permutations_tmm_normalized_cpm_no-cov.MEL_PBS.txt.gz --commands 40 commands.cpm.40.txt
while read c; do
     echo ${c} |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.cpm.no-cov.out" |\
     echo "submitted job"
done < commands.cpm.40.txt

###########################
# cpm of IFNg stimulated melanocytes (28 samples) to detect eQTL
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF4} --permute 1000 --normal --cov ${covF4} --out permutations_tmm_normalized_cpm.MEL_IFN.txt.gz --commands 40 commands.ifn.cpm.40.txt
while read c; do
     echo ${c} |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.ifn.cpm.out" |\
     echo "submitted job"
done < commands.ifn.cpm.40.txt

###########################
# fold-change to detect responseQTL
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF5} --permute 1000 --normal --cov ${covF5} --out permutations_FC_p0.05.MEL.txt.gz --commands 40 commands.FC.40.txt
while read c; do
     echo ${c} |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.FC.out" |\
     echo "submitted job"
done < commands.FC.40.txt

# combine all output together
cat permutations_tmm_normalized_log2cpm.MEL_PBS.txt.gz* | gzip -c > permutations_tmm_normalized_log2cpm_MEL_PBS.txt.gz
cat permutations_tmm_normalized_cpm.MEL_PBS.txt.gz* | gzip -c > permutations_tmm_normalized_cpm.MEL_PBS.txt.gz
cat permutations_log2FC_p0.05.MEL.txt.gz* | gzip -c > permutations_log2FC_p0.05.MEL.txt.gz
cat permutations_tmm_normalized_cpm_no-cov.MEL_PBS.txt.gz* | gzip -c > permutations_tmm_normalized_cpm_no-cov.MEL_PBS.txt.gz
rm permutation*.txt.gz.chr*

##############################
# filter genotype and dosage tables based on selected genes
# purpose: reduce memory required to load genotype table in R
cat permutations_*.txt | cut -f6 | sort | uniq > filtered.SNPs.txt
split -d -l 10000 filtered.SNPs.txt filtered.SNPs. 
f=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.genotype.txt #option2: filtered.dosage.txt
#f=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.dosage.txt
head -1 ${f} | sed 's/\[[0-9]*\]//g' | sed 's/:GT//g' | sed 's/\# //g' > temp.header
tail -n +2 ${f} > temp.body
split -d -l 100000 temp.body temp.split.
for f in temp.split.*;do
  for idx in 00 01 02 03 04;do
    grep -f filtered.SNPs.${idx} ${f} > grepped.${f}.${idx}
    echo "done with" ${f} ${idx}
    date
  done
done
cat temp.header grepped.* > fastQTL_selected.genotype.txt
#cat temp.header grepped.* > fastQTL_selected.dosage.txt
rm temp.* grepped.*


#############################
# 09/27/2021
# running fastQTL for selected genes: IRF4, IL15, IL15RA
# region + and - 10Kb
region_IRF4=chr6:381752-421443
region_IL15RA=chr10:5942372-5987611
region_IL2RA=chr10:6012071-6062309

# IRF4
fastQTL --vcf ${genotypeF} --bed ${phenotypeF3} --region chr6:300000-450000 --permute 1000 --out permutations.quantile.IRF4.log2FC.txt.gz --normal --cov ${covF3}
fastQTL --vcf ${genotypeF} --bed ${phenotypeF3} --region chr6:300000-450000 --out nominals.quantile.IRF4.log2FC.txt.gz --normal --cov ${covF3}
fastQTL --vcf ${genotypeF} --bed ${phenotypeF2} --region chr6:300000-450000 --permute 1000 --out permutations.quantile.IRF4.log2cpmPBS.txt.gz --normal --cov ${covF2}
fastQTL --vcf ${genotypeF} --bed ${phenotypeF2} --region chr6:300000-450000 --out nominals.quantile.IRF4.log2cpmPBS.txt.gz --normal --cov ${covF2}
# IL15RA
fastQTL --vcf ${genotypeF} --bed ${phenotypeF3} --region chr10:5942372-5987611 --permute 1000 --out permutations.quantile.IL15RA.log2FC.txt.gz --normal --cov ${covF3} 
fastQTL --vcf ${genotypeF} --bed ${phenotypeF3} --region chr10:5942372-5987611 --out nominals.quantile.IL15RA.log2FC.txt.gz --normal --cov ${covF3} 
# not log2cpm PBS data
# IL2RA
fastQTL --vcf ${genotypeF} --bed ${phenotypeF3} --region chr10:6000000-6100000 --permute 1000 --out permutations.quantile.IL2RA.log2FC.txt.gz --normal --cov ${covF3}
fastQTL --vcf ${genotypeF} --bed ${phenotypeF3} --region chr10:6000000-6100000 --out nominals.quantile.IL2RA.log2FC.txt.gz --normal --cov ${covF3}
# not log2cpm PBS data nor log2FC data


###############################
# 10/11/2021
# running fastQTL for skin_disease_GWAS_SNPs filtered genotype data
module load boost/1.75.0
module load gsl/2.6
module load g++/8.1.0
module load binutils/2.37
export PATH=${PATH}:/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/filtered_skin_disease_GWAS_SNPs.vcf.gz
phenotypeF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm.MEL_PBS.bed.gz
phenotypeF3=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log2FC_p0.05.MEL.bed.gz
covF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_MEL_PBS_covariates.txt
covF3=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log2FC_p0.05_MEL_covariates.txt

###########################
# TMM-normalized log2 CPM from Melanocyte-PBS samples
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF2} --permute 1000 --normal --cov ${covF2} --out permutations_GWASSNPs_tmm_normalized_log2cpm.MEL_PBS.txt.gz --commands 40 commands.40.txt
while read c; do
     echo "${c}" |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10000] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.GWASSNPs_tmm_normalized_log2cpm.MEL_PBS.out"
done < commands.40.txt
# make sure to rename files

###########################
# log2 fold change from Melanocyte samples IFNg vs PBS
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF3} --permute 1000 --normal --cov ${covF3} --out permutations_GWASSNPs_log2FC.MEL.txt.gz --commands 40 commands.40.txt
while read c; do
     echo "${c}" |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10000] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.GWASSNPs_tmm_normalized_log2cpm.MEL_PBS.out"
done < commands.40.txt
# make sure to rename files
cat permutations_GWASSNPs_tmm_normalized_log2cpm.MEL_PBS.txt.gz* | gzip -c > permutations_GWASSNPs_tmm_normalized_log2cpm_MEL_PBS.txt.gz
#####
# this did not result in anything.



###########################
# I need to prune the genotype file with MAF > 0.01
module load bcftools/1.9
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.vcf.gz
bcftools view -q 0.01:minor ${genotypeF} -Oz -o filtered.MAF0.01.vcf.gz
bcftools view -q 0.05:minor ${genotypeF} -Oz -o filtered.MAF0.05.vcf.gz
bcftools view -q 0.20:minor ${genotypeF} -Oz -o filtered.MAF0.20.vcf.gz

module load boost/1.75.0
module load gsl/2.6
module load g++/8.1.0
module load binutils/2.37
export PATH=${PATH}:/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/filtered.MAF0.05.vcf.gz
phenotypeF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm.MEL_PBS.bed.gz
covF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_MEL_PBS_covariates.txt

# TMM-normalized log2 CPM from Melanocyte-PBS samples
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF2} --permute 1000 --normal --cov ${covF2} --out permutations_MAF0.05_tmm_normalized_log2cpm.MEL_PBS.txt.gz --commands 40 commands.40.txt
while read c; do
     echo "${c}" |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10000] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.MAF0.05_tmm_normalized_log2cpm.MEL_PBS.out"
done < commands.40.txt
cat permutations_MAF0.05_tmm_normalized_log2cpm.MEL_PBS.txt.gz* | gzip -c > permutations_MAF0.05_tmm_normalized_log2cpm.MEL_PBS.txt.gz
rm permutations_MAF0.05_tmm_normalized_log2cpm.MEL_PBS.txt.gz.chr*
#### checking using fastQTL_melanocytes.Rmd, I realized this resulted in less power!
#### so fewer variants = less power.

##############################
# subset to only a few genes

module load condas/2018-05-11
source activate sshan_isoform
module load tabix/0.2.6
module load boost/1.75.0
module load gsl/2.6
module load g++/8.1.0
module load binutils/2.37

# make covariate table
ln -s /nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/tmm_normalized_log2cpm_top2500_variable.MEL_PBS.bed
bgzip tmm_normalized_log2cpm_top2500_variable.MEL_PBS.bed && tabix -p bed tmm_normalized_log2cpm_top2500_variable.MEL_PBS.bed.gz
source activate peer
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
Rscript run_PEER.R tmm_normalized_log2cpm_top2500_variable.MEL_PBS.bed.gz tmm_normalized_log2cpm_top2500_variable.MEL_PBS 15
cat genotype_dosage_topPCs.txt tmm_normalized_log2cpm_top2500_variable.MEL_PBS.PEER_covariates.txt > tmm_normalized_log2cpm_top2500_variable_MEL_PBS_covariates.txt # manually tweak the column name by adding id to the front.

genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.vcf.gz
phenotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_top2500_variable.MEL_PBS.bed.gz
#covF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_top2500_variable_MEL_PBS_covariates.txt
covF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_MEL_PBS_covariates.txt

# TMM-normalized log2 CPM from Melanocyte-PBS samples
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
export PATH=${PATH}:/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin
fastQTL --vcf ${genotypeF} --bed ${phenotypeF} --permute 1000 --normal --cov ${covF} --out permutations_tmm_normalized_log2cpm_top2500_variable_MEL_PBS.txt.gz --commands 40 commands.40.txt
while read c; do
     echo "${c}" |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10000] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.tmm_normalized_log2cpm_top2500_variable.MEL_PBS.out"
done < commands.40.txt
cat permutations_tmm_normalized_log2cpm_top2500_variable_MEL_PBS.txt.gz* > permutations_tmm_normalized_log2cpm_top2500_variable_MEL_PBS.txt
rm permutations_tmm_normalized_log2cpm_top2500_variable_MEL_PBS.txt.gz.chr*


#################################
##### 10/12/2021
### run fastQTL with different window sizes instead of 1Mb window default. Goal: see how that changes power
### Options:
# 1. 500Kb. prefix=tmm_normalized_log2cpm_MEL_PBS_500Kb_window
# 2. 100Kb. prefix=tmm_normalized_log2cpm_MEL_PBS_100Kb_window
# 3. 2Mb.   prefix=tmm_normalized_log2cpm_MEL_PBS_2Mb_window
module load condas/2018-05-11
source activate sshan_isoform
module load gsl/2.6
module load g++/8.1.0
module load binutils/2.37

genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.vcf.gz
phenotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm.MEL_PBS.bed.gz
covF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_MEL_PBS_covariates.txt
prefix=tmm_normalized_log2cpm_MEL_PBS_100Kb_window
size=1e5

# TMM-normalized log2 CPM from Melanocyte-PBS samples
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
export PATH=${PATH}:/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin
fastQTL --vcf ${genotypeF} --bed ${phenotypeF} --permute 1000 --normal --window ${size} --cov ${covF} --out temp.permutations_${prefix}.txt --commands 40 commands.40.txt
while read c; do
     echo "${c}" |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10000] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutations_${prefix}.out"
done < commands.40.txt
# once jobs are finished, run this:
cat temp.permutations_${prefix}.txt* > permutations_${prefix}.txt
rm temp.permutations_${prefix}.txt.chr*



#################################
##### 10/21/2021
### run fastQTL with a subset of SNPs (GWAS SNPs for skin diseases, and variants in high LD (r2>0.6) in 1Mb distance)
module load condas/2018-05-11
source activate sshan_isoform
module load gsl/2.6
module load g++/8.1.0
module load binutils/2.37
module load bcftools/1.9

cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/GWAS_SNP_with_LD
# genotypeF generated by /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/log.sh 

genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/GWAS_SNP_with_LD/filtered_skin_disease_GWAS_SNPs.vcf.gz
phenotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm.MEL_PBS.bed.gz
covF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_MEL_PBS_covariates.txt
prefix=tmm_normalized_log2cpm_MEL_PBS_SNPs_with_LD

export PATH=${PATH}:/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin
fastQTL --vcf ${genotypeF} --bed ${phenotypeF} --include-sites ${subset} --cov ${covF} --permute 1000 --normal --out temp.permutations_${prefix}.txt --commands 40 commands.40.txt
while read c; do
     echo "${c}" |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10000] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutations_${prefix}.out"
done < commands.40.txt
# once jobs are finished, run this:
cat temp.permutations_${prefix}.txt* > permutations_${prefix}.txt
cat temp | grep -v "NA" > permutations_${prefix}.txt
cat temp | grep "NA" | cut -f1 > genes_without_variants.txt
rm temp
rm temp.permutations_${prefix}.txt.chr*








