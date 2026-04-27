#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=2040]
#BSUB -q long
#BSUB -W 08:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

# set-up working directory
dir=/pi/manuel.garber-umw/human/skin/eQTLs/edQTL

# step 1. compute_editing_level.sh
# 01/22/2023 the coverage from genotyping bam files is too low (around 0-2 reads per site). 
# 01/24/2023 trying to see high coverage from MEL RNAseq PBS bam file. 4,768 editing sites with coverage > 10. 2,741 editing sites with coverage > 20. 113 editing sites with coverage > 20 and editing level > 0. 136 editing sites with coverage > 10 and editing level > 0
# 01/24/2023 try with merging PBS and IFN bam files. 6,594 editing sites with coverage > 10. 3,592 editing sites with coverage > 20. 149 editing sites with coverage > 20 adn editing level > 0, 195 editing sites with coverage > 10 and editing level > 0.
# 01/24/2023 try with merging all bam files from all 3 celltypes. Since the genotyping result for each celltype should be the same. 31,113 editing sites with coverage > 10. 17,421 editing sites with coverage > 20. 1,157 sites have coverage > 10 and editing level > 0. 874 sites have coverage > 20 and editing level > 0. 
# ^ I will go with the option above, merging all RNAseq bam files available for each donor.
cd ${dir}/output
while read c; do
  echo ${c} | bsub -J edQTL -W 5:00 -n 1 -R "span[hosts=1]" -R rusage[mem=4080] -q long -e "${dir}/log/step1.%J%I.err" -o "${dir}/log/step1.%J%I.out"
done < ${dir}/scripts/commands_compute_editing_level.txt

# step 2. combine individual editing data into a matrix
cd ${dir}/output
conda activate fastQTL
perl ${dir}/scripts/sharedsamples_sites_matrix_FastQTL_v8.pl > foreskin.edMat.10cov.20samps.txt
head -1 foreskin.edMat.10cov.20samps.txt > foreskin.edMat.10cov.20samps.noXYM.txt
cat foreskin.edMat.10cov.20samps.txt | grep -E "chr[0-9]{1,2}" >> foreskin.edMat.10cov.20samps.noXYM.txt
conda deactivate

# step 3. convert editing level matrices to format recognized by edQTL
conda activate fastQTL
python ${dir}/scripts/prepare_phenotype_table_for_QTLtools.V8.py --pcs 15 foreskin.edMat.10cov.20samps.noXYM.txt
sh foreskin.edMat.10cov.20samps.noXYM.txt_prepare.sh
echo -e "#chr"'\t'"start"'\t'"end"'\t'"name" > header1
head -1 foreskin.edMat.10cov.20samps.noXYM.txt | tr ' ' '\t' | cut -f2- | paste header1 - > header.txt
zcat foreskin.edMat.10cov.20samps.noXYM.txt.qqnorm_chr*gz | grep -v "start" | sed 's/nan/NA/g' | sort -k 1,1 -k2n,2 | awk '{print "chr"$1,$2,$3,"chr"$1"_"$3,$0}' OFS="\t" | cut -f 1,2,3,4,9- | sed 's/^chr#/#/g' > body.txt
cat header.txt body.txt | bgzip -c > foreskin.edMat.10cov.20samps.noXYM.qqnorm.bed.gz
rm header1 header.txt body.txt
tabix -f -p bed foreskin.edMat.10cov.20samps.noXYM.qqnorm.bed.gz
conda deactivate

# step 4. Obtain PEER factors from the combined editing level matrices
conda activate peer # peer is its own conda env because it runs on older version of R: 3.4.1
Rscript ${dir}/scripts/run_PEER.R \
    foreskin.edMat.10cov.20samps.noXYM.qqnorm.bed.gz \
    foreskin.edMat.10cov.20samps.noXYM.qqnorm \
    10
conda deactivate

# step 5. Combine genotype PCs with phenotype PEER factors
conda activate fastQTL
genotypePCs=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/pca/pca_eigenvec_table.txt
Rscript ${dir}/scripts/combine_covariates.R \
    foreskin.edMat.10cov.20samps.noXYM.qqnorm.PEER_covariates.txt \
    ${genotypePCs} \
    ${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.combined_covariates.txt
bgzip -c foreskin.edMat.10cov.20samps.noXYM.qqnorm.combined_covariates.txt > foreskin.edMat.10cov.20samps.noXYM.qqnorm.combined_covariates.txt.gz
conda deactivate

# step 6. Running FastQTL: nominal step
conda activate fastQTL
export PATH=$PATH:/pi/manuel.garber-umw/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin
vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
phenotypeF=${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.bed.gz
covF=${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.combined_covariates.txt.gz
jobname=fastQTL
fastQTL --vcf ${vcfF} --bed ${phenotypeF} --cov ${covF} --out foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal --permute 1000 --window 1e5 --commands 40 commands.40.txt
while read c; do
	echo "
	  export PATH=$PATH:/pi/manuel.garber-umw/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin;
	  dir=/pi/manuel.garber-umw/human/skin/eQTLs/edQTL;
	  vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz;
	  phenotypeF=${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.bed.gz;
	  covF=${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.combined_covariates.txt.gz;
	  ${c}" |\
	bsub -J ${jobname} -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10000] -q long -o "${dir}/log/nominal.%J%I.out" -e "${dir}/log/nominal.%J%I.err" |\
	echo "submitted job"
done < commands.40.txt	  
while [[ $(bjobs | grep ${jobname} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${jobname} | wc -l) "jobs remaining"; date;sleep 5; done
cat foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal.chr* > ${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal
rm foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal.chr*
conda deactivate

# clean-up
rm foreskin.edMat.10cov.20samps.noXYM.txt.phen_chr*
rm foreskin.edMat.10cov.20samps.noXYM.txt.qqnorm_chr*gz*
rm *.out

# step 7. Filter out significant results to plot in R
cd ${dir}/output
cat ${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal | awk '$11<0.05' > foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal.sig
cat foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal.sig | cut -d' ' -f6 | sort | uniq > temp.snps
cat foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal | cut -d' ' -f6 | sort | uniq > temp.snps
vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
module load bcftools/1.16
bcftools view -h ${vcfF} > temp.header.txt
bcftools view --include ID==@temp.snps ${vcfF} -Oz -o temp.snps.vcf.gz
bcftools view -H temp.snps.vcf.gz | cat temp.header.txt - > temp.snps.vcf.with.header
bcftools view temp.snps.vcf.with.header -Oz -o snps.vcf.gz
rm temp.header.txt temp.snps.vcf.with.header

# obtain the genotype and info from vcf file
genotype_table=${dir}/output/snps.vcf.gz
bcftools query -f "%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT{0}[\t%GT]\n" -H ${genotype_table} | head -1 > header
cat header | tr '\t' '\n' | cut -d']' -f2 | tr '\n' '\t'  | sed '$s/\t$/\n/' > header2
cat header2 | sed 's/POS\tPOS/START\tEND/g'  > header3
bcftools query -f "%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT{0}[\t%GT]\n" ${genotype_table} -o temp.filtered.genotype.vcf
cat temp.filtered.genotype.vcf | awk '{OFS="\t"}{print $1,$2,$2+1,$4,$5,$6}' > temp1
cat temp.filtered.genotype.vcf | cut -f7- > temp2
paste temp1 temp2 > temp3
cat temp3 > temp.filtered.genotype.bed
cat header3 > temp.header
cat temp.header temp.filtered.genotype.bed > genotype.bed
rm header header2 header3 temp.filtered.genotype.vcf temp1 temp2 temp3 temp.header temp.filtered.genotype.bed

# overlap edSites with genes
module load bedops/2.4.41
module load bedtools/2.30.0
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks 
cat Homo_sapiens.GRCh38.105.gtf | gtf2bed - > Homo_sapiens.GRCh38.105.gtf.bed
cat Homo_sapiens.GRCh38.105.gtf.bed | awk '$8=="gene"' | awk '{print "chr"$0}' > Homo_sapiens.GRCh38.105.genes.gtf.bed
cd ${dir}/output
cat ${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal | cut -d' ' -f1 | awk -F"_" '{OFS="\t"}{print $1,$2-1,$2,$1"_"$2}' | sort -k 1,1 -k2,2n > temp.edSites
cat /pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/expressedGenes.txt \
    /pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/expressedGenes.txt \
    /pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/expressedGenes.txt |\
    sort | uniq | awk '{print "\x22"$0"\x22"}' >  expressedGenes.txt
grep -f expressedGenes.txt /pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks/Homo_sapiens.GRCh38.105.genes.gtf.bed > GRCh38.105.expressedGenes.gtf.bed
edSites=${dir}/output/temp.edSites
expressedGenes=${dir}/output/GRCh38.105.expressedGenes.gtf.bed
bedtools closest -d -a ${edSites} -b ${expressedGenes} > temp1
cat temp1 | awk '{print $4}' > temp2.edSite
cat temp1 | cut -f14 | sed 's/.*gene_name \"//g' | cut -d'"' -f1 > temp2.gene
cat temp1 | cut -f15 > temp2.dist
paste temp2.edSite temp2.gene temp2.dist > edsites_closest_expressed_genes.txt
rm temp1 temp2.*


# plot in R
conda activate fastQTL
zcat ${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.bed.gz > ${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.bed
edSiteCountF=${dir}/output/foreskin.edMat.10cov.15samps.txt
edSiteLevelF=${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.bed
genotypeBedF=${dir}/output/genotype.bed
fastQTLResF=${dir}/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal.sig
Rscript ${dir}/scripts/paired_plot.R ${dir} ${fastQTLResF} ${edSiteCountF} ${edSiteLevelF} ${genotypeBedF}
















