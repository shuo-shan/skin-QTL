#!/bin/bash
#BSUB -n 9
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=45000]
#BSUB -q long
#BSUB -W 121:00
### script to interface Gencove to retrieve, upload, and process genotyping data

### crystal shan 06/2021
### bsub -Is -q interactive -W 8:00 -n1 -R rusage[mem=450000] -R "span[hosts=1]" /bin/bash

############################
#### download gencove data in hg38
module load condas/2018-05-11
source activate sshan_isoform
#pip install gencove
##To get that run through, you'd need to do the following:
##- Download the fastq files from Gencove
##- Reupload them under a new data analysis configuration
##- Run the pipeline
#### check and restore archived samples
#gencove projects list
#gencove projects list-samples 481d8927-d82c-4865-b00b-d530f346041c
#gencove projects restore-samples 481d8927-d82c-4865-b00b-d530f346041c
#### download fastq
#gencove download . --project-id 481d8927-d82c-4865-b00b-d530f346041c --file-types fastq-r1,fastq-r2,metadata
#### check md5sum and backup to Amazon
#cd /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/fastq
#for i in `find . -type f`; do md5sum $i > $i.md5sum; done
#/project/umw_biocore/bin/amazonBackup.bash fastq s3://biocorebackup/garberlab/human/skin/eQTL/genotyping
#### upload fastq
#gencove upload . gncv://batch-1/
#### download hg38-processed data
#cd /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf 
#gencove download . --project-id bd841628-fcc2-487a-8460-f5428237f0c9 --file-types impute-vcf
#
############################
###### process vcf files
module load bcftools/1.9
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/organized
# set variants with LOWCONF (max(GP)<0.90) to have missing genotype (GT=./.)
#ls $dir | grep '.vcf.gz' | grep -v 'csi' | tr '\t' '\n' | grep -v 'md5sum' > $dir/filelst.txt
#outdir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/organized_lowconfGTmod
#cd ${outdir}
#while read donor;do
#  echo "working on" ${donor}
#  date
#  inf=${dir}/${donor}.vcf.gz
#  echo $inf
#  outf=${outdir}/${donor}.vcf.gz
#  echo ${donor} > this.txt
#  bcftools +setGT --threads 8 ${inf} -- -t q -i 'max(GP)<0.90' -n "./." > this.vcf
#  bcftools reheader -s this.txt --threads 8 this.vcf -o this2.vcf
#  bcftools view this2.vcf -v snps --threads 8 -Oz -o ${outf}
#  rm this.txt this.vcf this2.vcf
#  echo "done with" ${donor}
#  date
#done < $dir/filelst.txt
## index all .vcf.gz files
dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/organized_lowconfGTmod
cd $dir
for f in *.vcf.gz;do
  bcftools index --threads 8 ${f}
  date +"[%d-%m-%y] %T"
  echo "indexing done for " ${f}
done
# merge all .vcf.gz files
date
echo "begin to merge all donor files"
# manually filter out filelst entries that don't have RNAseq data.
dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/organized_lowconfGTmod
ls $dir | grep '.vcf.gz' | grep -v 'csi' | grep -v 'CB' | grep -v 'VB' > $dir/filelst.txt
outdir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged
cd ${dir}
bcftools merge -l $dir/filelst.txt -Oz -o $outdir/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.vcf.gz --threads 8
date
echo "merging done. :)"
#rm $dir/filelst.txt
## filter merged .vcf.gz file
#cd $dir
## project id on Gencove
#f=bd841628-fcc2-487a-8460-f5428237f0c9
## filter by minimum alt allele read count
bcftools view $outdir/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.vcf.gz --min-ac 3 --threads 8 -Oz -o $outdir/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.vcf.gz
#bcftools view bd841628-fcc2-487a-8460-f5428237f0c9.merged.vcf.gz --min-ac 1 -Oz -o ${f}.merged.filtered.1.vcf.gz
echo "filter 1 done! :)"
date +"[%d-%m-%y] %T"
## add AF tag
bcftools +fill-tags $outdir/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.vcf.gz -Oz -o $outdir/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz -- -t AF
#bcftools +fill-tags ${f}.merged.filtered.1.vcf.gz -Oz -o ${f}.merged.filtered.1.AFtagged.vcf.gz -- -t AF
echo "AF tag done! :)"
date +"[%d-%m-%y] %T"
## filter for QC
#bcftools view ${f}.merged.filtered.1.AFtagged.vcf.gz -f .,PASS -Oz -o ${f}.merged.filtered.2.AFtagged.vcf.gz
#echo "filter 2 done! :)"
#date +"[%d-%m-%y] %T"
#bcftools index ${f}.merged.filtered.1.vcf.gz
#bcftools index ${f}.merged.filtered.1.AFtagged.vcf.gz
#bcftools index ${f}.merged.filtered.2.AFtagged.vcf.gz
bcftools index $outdir/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
echo "index done! :)"
date +"[%d-%m-%y] %T"
# indexing with tabix to feed into fastQTL
module load condas/2018-05-11
source activate sshan_isoform
conda install -c bioconda tabix
module load tabix/0.2.6
tabix -p vcf $outdir/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
#############################
## phasing 
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged
#f=bd841628-fcc2-487a-8460-f5428237f0c9
#inf=${f}.merged.filtered.2.AFtagged.vcf.gz
#cd $dir
#module load beagle/5.0
#module load java/1.8.0_171
#date +"[%d-%m-%y] %T"
#java -jar /share/pkg/beagle/5.0/beagle.16May19.351.jar gt=${inf} out=${dir}/${f}.phased iterations=5 nthreads=8
#date +"[%d-%m-%y] %T"
#echo "phasing done :)"
#bcftools index ${dir}/${f}.phased.vcf.gz
############################
## renaming sample id to name
#zcat bd841628-fcc2-487a-8460-f5428237f0c9.phased.vcf.gz.vcf.gz | sed 's/1912f712-8173-4b3d-a435-b8823226d0ec/CB032/g; s/381a63a5-fb6e-4fac-a59a-bc16d50f5209/CB043/g; s/9b485baa-f48a-4ded-a0a4-cadf877174c3/CB045/g; s/4482070f-8aea-41bc-aae6-d7314fc19eb2/F22/g; s/ef53c18a-f2f5-4fd5-b958-3be322c48b2f/F23/g; s/7091acee-30cd-4778-aeae-950dbc032cd4/F24/g; s/6d0b04d1-84e3-4979-8b8b-b93792da776f/F25/g; s/149430e7-6d2d-4e50-baba-3e8f044ce2e3/F27/g; s/e64cf0cc-e528-4ee7-8876-6e1b471f7c3f/F28/g; s/5b67f1d5-29b4-448f-b9a3-d23a2d12cd93/F30/g; s/d8380c80-b025-4d9b-9d06-67ee06eac4ad/F31/g; s/231efd05-64a1-4243-b292-0e454cb93de1/F32/g; s/0803da36-0f14-4df2-b68f-be1424566f2e/F33/g; s/5df7279b-bef9-44b3-a1c1-87b3318b79d4/F34/g; s/89ba837c-8125-4937-a1f4-04491da9b394/F35/g; s/93cfd69a-1288-4562-8894-5e8d3afca153/F36/g; s/61f190c4-cec6-45d9-8a70-24b3aa3a7902/F37/g; s/bb55dd6b-9f98-414f-a72c-856e3faa5243/F38/g; s/442a0e66-c794-46c1-aae5-d5d474e344c5/F39/g; s/6e085efe-ff90-48df-bc6f-08e33af14a5c/F40/g; s/e7b548bc-bed5-4a30-b4d9-0151966e6b9d/F41/g; s/30ace0f7-cb78-4da0-ada5-beff07ab1b9c/F42/g; s/ccc4e5cc-bb4f-4d9a-b5f2-853c6fce5623/F44/g; s/f4a987ea-ca08-4b38-8668-48e5efe45076/F45/g; s/f0b571af-c87e-4611-8288-bfcf4f9e9133/F46/g; s/a6ba9a8c-f89d-4a6e-b85d-05f0d0267d5d/F47/g; s/7d7cc250-d397-4ef4-97d8-94e7a2a85d8d/F48/g; s/b0690c72-9cf6-465f-b5dc-e0b94fe8c8f3/F49/g; s/edbb9163-65d7-4fc8-9d14-627d0987f087/F50/g; s/57948dd0-a324-4477-a917-c33f9c4e2daf/F51/g; s/383fe24c-a10e-424e-abcf-29d1e2a61c5d/F52/g; s/3ae0881d-e215-48be-8a5d-0f81736e715d/F53/g; s/8f482258-bffd-40ca-be3a-487644f0c52f/F55/g; s/29026b76-493d-4f87-9e3b-800e311c19ed/F56/g; s/b0a0708d-c71d-4f0b-8404-00bd03ee6b63/F57/g;  s/58993494-a077-4b70-866f-edcd06a79012/F58/g; s/3e47eefb-c67f-4814-bc33-3d3e616daeef/F59/g; s/49d91fb6-add6-497a-8c51-c95e5df25f0a/F60/g; s/839f2dfa-a0d9-46bf-9588-21c78c55b2d0/F61/g; s/3e5e3f37-fe51-4fb7-a356-80d803c3db9b/F62/g; s/413a4cbe-af1f-4cbb-a794-14c88b30da31/F63/g; s/6e008d11-42b2-4d3f-9bf3-ff78f1524cfa/VB126/g; s/9f5d36dd-208e-4a7f-b214-97d53020d6bb/VB150/g; s/9bb4e939-c6d5-4cae-854e-c5a184dedcc1/VB151/g; s/b3b92ecf-743f-4a02-b66c-f11f1fac7f2e/VB159/g; s/dec895e5-bf8b-4fce-8804-12fcd22bce5d/VB163/g; s/8f9db5b0-55fd-4236-a9b0-6f355178f9c9/VB172/g; s/9b9a8181-3236-40bc-a851-b1cfe456d2e2/VB173/g' > bd841628-fcc2-487a-8460-f5428237f0c9.phased.renamed.vcf
#bcftools view ${f}.phased.renamed.vcf -Oz -o ${f}.phased.renamed.vcf.gz
#
############################
###### process vcf files
#module load bcftools/1.9
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged
## filter merged .vcf.gz file
#cd $dir
## project id on Gencove
#f=bd841628-fcc2-487a-8460-f5428237f0c9
## filter by minimum alt allele read count
#date +"[%d-%m-%y] %T"
#bcftools view ${f}.merged.filtered.2.AFtagged.vcf.gz -H | cut -f8 | cut -d";" -f1 | cut -d"=" -f2 | uniq -c > AC.summary.txt
#echo "done summarizing allele counts"
#date +"[%d-%m-%y] %T"
###########################
###### process vcf file into matrix that we can cluster
#module load condas/2018-05-11
#source activate sshan_isoform
#RscriptPATH=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38
#Rscript ${RscriptPATH}/vcf_processing.r 
#
###########################
###### overlap phased filtered variants with chromatin accessible regions (ATACseq summits from 9bp shifted regions, extended)
#atacpeaks=/nl/umw_manuel_garber/kensei/skin_human/report4846/output/report4853/filteredbam/output_with_filter/report4853/atac/merged_peaks_slop200bp_concactnated.bed
#bcftools view -H bd841628-fcc2-487a-8460-f5428237f0c9.phased.vcf.gz -R ${atacpeaks} -Oz -o bd841628-fcc2-487a-8460-f5428237f0c9.phased.chromAcc.vcf.gz
#



#############################
## find merged vcf file for the 3 donors
#bcftools view -s 6d0b04d1-84e3-4979-8b8b-b93792da776f,b0690c72-9cf6-465f-b5dc-e0b94fe8c8f3,8f482258-bffd-40ca-be3a-487644f0c52f -Oz -o temp.vcf.gz bd841628-fcc2-487a-8460-f5428237f0c9.merged.filtered.2.AFtagged.vcf.gz
#bcftools view --min-ac 1 -Oz -o bd841628-fcc2-487a-8460-f5428237f0c9.filtered.2.3donors.vcf.gz temp.vcf.gz
#
## phasing 
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged
#f=bd841628-fcc2-487a-8460-f5428237f0c9
#inf=${f}.filtered.2.3donors.vcf.gz
#cd $dir
#module load beagle/5.0
#module load java/1.8.0_171
#date +"[%d-%m-%y] %T"
#java -jar /share/pkg/beagle/5.0/beagle.16May19.351.jar gt=${inf} out=${dir}/${f}.filtered.2.3donors.phased iterations=10 nthreads=8
#date +"[%d-%m-%y] %T"
#echo "phasing done :)"
#bcftools index ${dir}/${f}.filtered.2.3donors.phased.vcf.gz
#


#############################
# find merged vcf file for 16 donors with RNAseq MEL data
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged
#cd $dir
#module load bcftools/1.9
#bcftools view -s 149430e7-6d2d-4e50-baba-3e8f044ce2e3,5b67f1d5-29b4-448f-b9a3-d23a2d12cd93,d8380c80-b025-4d9b-9d06-67ee06eac4ad,5df7279b-bef9-44b3-a1c1-87b3318b79d4,e7b548bc-bed5-4a30-b4d9-0151966e6b9d,f4a987ea-ca08-4b38-8668-48e5efe45076,f0b571af-c87e-4611-8288-bfcf4f9e9133,b0690c72-9cf6-465f-b5dc-e0b94fe8c8f3,edbb9163-65d7-4fc8-9d14-627d0987f087,57948dd0-a324-4477-a917-c33f9c4e2daf,383fe24c-a10e-424e-abcf-29d1e2a61c5d,8f482258-bffd-40ca-be3a-487644f0c52f,29026b76-493d-4f87-9e3b-800e311c19ed,6d0b04d1-84e3-4979-8b8b-b93792da776f,3e5e3f37-fe51-4fb7-a356-80d803c3db9b,413a4cbe-af1f-4cbb-a794-14c88b30da31 --min-ac 1 -Oz -o bd841628-fcc2-487a-8460-f5428237f0c9.merged.filtered.2.16donors.vcf.gz bd841628-fcc2-487a-8460-f5428237f0c9.merged.filtered.2.AFtagged.vcf.gz 
#
## phasing
#f=bd841628-fcc2-487a-8460-f5428237f0c9
#inf=${f}.merged.filtered.2.16donors.vcf.gz
#cd $dir
#module load beagle/5.0
#module load java/1.8.0_171
#date +"[%d-%m-%y] %T"
#java -jar /share/pkg/beagle/5.0/beagle.16May19.351.jar gt=${inf} out=${dir}/${f}.merged.filtered.2.16donors.phased iterations=10 nthreads=8
#date +"[%d-%m-%y] %T"
#echo "phasing done :)"
#
##############################
### discard sites with more than 0.8 linkeage disequilibrium
#module load bcftools/1.9
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged
#f=bd841628-fcc2-487a-8460-f5428237f0c9
#inf=${dir}/${f}.merged.filtered.2.16donors.vcf.gz
#bcftools +prune -l 0.8 -w 50 ${inf} -Oz -o ${dir}/${f}.16donors.filtered2.pruned.vcf.gz
#

#########################
#### process vcf file into matrix that we can cluster
#module load condas/2018-05-11
#source activate sshan_isoform
#RscriptPATH=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38
#Rscript ${RscriptPATH}/vcf_processing.r 
#
###########################
## 07/19/2021
#module load bcftools/1.9
#cd /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/organized
#for f in *.vcf.gz;do
#  bcftools view -H ${f} | grep -m 1 -w "rs706779" >> /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/hg38.rs706779.grep.out
#  echo "done with "${f}
#done
#
##########################
#cd /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping 
#bcftools view -H filt_renamed_impute-vcf-merged_with_AF.vcf.gz | grep -m 1 -w "rs706779" >> /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/hg19.filtvcf.rs706779.grep.out







