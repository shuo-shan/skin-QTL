#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 250:00
#BSUB -e "./pipeline.%J%I.err"
#BSUB -o "./pipeline.%J%I.out"

# written by Crystal Shan 06/2022, modified 05/2024
# goal: take reQTLs and PBS and IFN eQTLs (lose cutoff pval<0.01), look for reQTL via beta-comparison model#8

# set-up
module load bedtools/2.30.0
celltype=KRT
jobname=${celltype}
dir=/pi/manuel.garber-umw/human/skin/eQTLs/case_studies/DHX58/rs8081327/model_betacomparison
cd ${dir}
mkdir plot log

# copy the covariate files, expressed gene files from basic DREG pipeline
olddir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/KRT_minimal
cp ${olddir}/expressed_genes_tss* ${dir}
cp ${olddir}/snps_near_expressed_genes.bed ${dir}

# get the eligible SNP:gene pair (400,149)
cd /pi/manuel.garber-umw/human/skin/eQTLs/case_studies/DHX58/rs8081327
cat KRT_DESeq2_results_IFNvsPBS.txt | tr ' ' '\t' | awk '{if ($8<0.05) print $0}' | awk 'NR==1{print $1,$2,$3,$4,$5,$6,$7}NR>1{if ($4<-1.5 || $4>1.5) print $2,$3,$4,$5,$6,$7,$8}' > KRT_DESeq2_results_IFNvsPBS_sigDEgenes_padj0.05absLog2FC1.5.txt
cat KRT_DESeq2_results_IFNvsPBS_sigDEgenes_padj0.05absLog2FC1.5.txt | awk 'NR>1{print $1}' | sort | uniq > KRT_DEgenes.txt
cd model_betacomparison/
cat ../KRT_DEgenes.txt | awk '{OFS=FS="\t"}{print "rs8081327",$0}' > SNP_gene_pairs.txt

# for each SNP:gene pair, run the beta-comparison model with 10,000 permutations
mkdir ${dir}/log/build_models
mkdir ${dir}/modelingResult
cat ${dir}/SNP_gene_pairs.txt | awk -v dir=${dir} '{print "bash "dir"/pipeline_betacomparison_KRT_perSNPGenePair.sh "$0" "dir" masteroutput"}' > ${dir}/commands.txt
sleep 3
bash /pi/manuel.garber-umw/sshan/scripts/function_collapse_commands.sh ${dir} ${dir}/commands.txt ${dir}/commands.joined.txt 10
sleep 5
while read c; do
        echo ${c} | bsub -W 72:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=800] -q long -e "${dir}/log/build_models_%J%I.err" -o "${dir}/log/build_models_%J%I.out"
done < ${dir}/commands.joined.txt

while [[ $(bjobs | grep ${jobname} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${jobname} | wc -l) "jobs remaining"; date;sleep 60; done



