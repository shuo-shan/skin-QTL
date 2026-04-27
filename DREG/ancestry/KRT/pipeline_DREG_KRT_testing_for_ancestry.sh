#!/bin/bash
# this script is run my DREG model to test whether ancestry is a feature in any snp-gene pairs, especially significant QTL-Genes.
#BSUB -n 1
#BSUB -J pipeline
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 250:00
#BSUB -e "./pipeline.%J%I.err"
#BSUB -o "./pipeline.%J%I.out"

# written by Crystal Shan 06/2022, modified 07/2023
# goal: take SNPs that are nearby any gene, look for Differentially Responsive (DR) genes by the SNP genotype
# set-up
module load bedtools/2.30.0
celltype=KRT
jobname=${celltype}
dir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/ancestry/KRT
cd ${dir}
#mkdir plot log
#
#
####### 1. compile a bed file for all SNPs that are within 500Kbp range of an expressed autosomal gene's TSS.
## compile a bed file of 1Mbp windows around expressed gene's TSS
#datadir=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis
#tss=/pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks/Ensembl_GRCh38.105_transcription_start_sites.bed
#genome_sizes=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
#bash /pi/manuel.garber-umw/sshan/scripts/function_rapid_fgrep.sh ${datadir}/expressedGenes.txt ${tss} expressed_genes_tss.bed haha 1000
#cat expressed_genes_tss.bed | grep -E '\bchr[0-9]{1,2}\b' |\
#               bedtools sort -i stdin | bedtools slop -i stdin -b 500000 -g ${genome_sizes} |\
#               bedtools merge -i stdin | awk '{OFS="\t"}{print $0}' > expressed_genes_tss_flanking_merged.bed
#
## Filter SNPs to keep those within 500Kbp range from the TSS of an expressed gene
## note: the 1,924,039 genotype SNPs are already filtered to include SNPs that has 3 genotypes and at least 3 donors for each genotype
##       bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/vcf_filter_then_to_bed.sh
#head -1 /pi/manuel.garber-umw/human/skin/eQTLs/DREG/master_filtered_genotype.bed > temp.header
#bedtools intersect -a /pi/manuel.garber-umw/human/skin/eQTLs/DREG/master_filtered_genotype.bed \
#                  -b expressed_genes_tss_flanking_merged.bed |\
#                  cat temp.header - > snps_near_expressed_genes.bed
#
## clean up
#rm temp.header temp_expressed_genes.txt; rm -r rapid_fgrep_temp
## ----> these steps result in  1,454,078 snps that are near a gene's TSS, has at least 3 genotypes and at least 3 donors for each genotype
## these snps are not pruned because I want to prioritize SNPs within regions of high LD by incorporating functional genomic data later on.
#
## Optional: Further pruning SNP based on autosome or high LD.
##### option1: no further filtering, move on with the full list of SNPs
#cat snps_near_expressed_genes.bed | awk 'NR>1{print $0}' | cut -f4 > snps.txt
##### option2: keep autosome SNPs in low LD.
##source activate skineqtl # this runs +prune plugin
### pick all SNPs from autosomes (  524406 snps --> 509,360 snps)
##cat /pi/manuel.garber-umw/human/skin/eQTLs/DREG/${celltype}_filtered.genotype.bed | awk '{if ($1!="chrX") print $4}' | sort | uniq > snps_autosome.txt
### prune SNPs (509,360 snps --> 54,902 snps)
##bcftools view --include ID==@snps_autosome.txt ${genotype_table} -Oz -o autosome_snps.vcf.gz
###source activate vcflib_1.0.3 # this runs +prune plugin
##bcftools +prune -l 0.6 -w 1000 ${genotype_table} -Oz -o pruned_snps.vcf.gz
##bcftools view --include ID==@snps_autosome.txt pruned_snps.vcf.gz -Oz -o pruned_autosome_snps.vcf.gz
##bcftools index pruned_autosome_snps.vcf.gz
##bcftools view -H pruned_autosome_snps.vcf.gz | cut -f3 > snps.txt
#
####### 2. Compile log2FC covariate matrix based on the desired number of genotypePCs and latent phenotype variables.
#for n_PCs in 10;do
#    for n_latentVar in 10;do
#      mkdir ${dir}/modelingResult_a${n_PCs}b${n_latentVar}
#      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} CPM_PBS
#      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} CPM_IFN
#      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} Log2FCwithDummy10
#      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} rankNormCPM_PBS
#      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} rankNormCPM_IFN
#      bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/compile_covariates_KRT.sh ${n_PCs} ${n_latentVar} rankNormLog2FCwithDummy10
#    done
#done

######### 3. For every SNP of interest, find the nearby expressed genes, and construct the model
#mkdir ${dir}/log/build_models
#shuf snps.txt > shuffled_snps.txt
#mkdir split_snps; cd split_snps; split -l 50000 ${dir}/shuffled_snps.txt splitSnps_; rm ${dir}/shuffled_snps.txt
n_PCs=10
n_latentVar=10
cd ${dir}/split_snps
for this_file in now_splitSnps_*;do
        rm ${dir}/commands.txt ${dir}/shuffled_commands.txt ${dir}/commands.joined.txt
        while read snp;do
            echo "bash ${dir}/pipeline_DREG_KRT_perSNP_testing_for_ancestry.sh ${snp} a${n_PCs}b${n_latentVar} ${dir} masteroutput_a${n_PCs}b${n_latentVar}_${this_file}" >> ${dir}/commands.txt
        done < ${dir}/split_snps/${this_file}
        sleep 5
        shuf ${dir}/commands.txt > ${dir}/shuffled_commands.txt # randomly shuffle the rows (snps) so snps in gene-dense regions don't get stuck in the same joined commands. it significantly prolongs job time.
        bash /pi/manuel.garber-umw/sshan/scripts/function_collapse_commands.sh ${dir} ${dir}/shuffled_commands.txt ${dir}/commands.joined.txt 55
        sleep 5
        while read c; do
          echo ${c} | bsub -W 05:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=400] -q long -e "${dir}/log/build_models_%J%I.err" -o "${dir}/log/build_models_%J%I.out"
        done < ${dir}/commands.joined.txt
        sleep 5

        while [[ $(bjobs | grep ${jobname} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${jobname} | wc -l) "jobs remaining"; date;sleep 120; done

        # clean-up
        rm ${dir}/log/build_models_* 
        rm ${dir}/masteroutput_a${n_PCs}b${n_latentVar}_${this_file}.lock
done


