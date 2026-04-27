#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 250:00
#BSUB -e "./pipeline.%J%I.err"
#BSUB -o "./pipeline.%J%I.out"

# written by Crystal Shan 06/2022, modified 08/21/2024
# set-up
module load bedtools/2.30.0
celltype=KRT
jobname=${celltype}
dir=/pi/manuel.garber-umw/human/skin/eQTLs/case_studies/DHX58/rs8081327/model_log2FClinear
cd ${dir}
mkdir plot log

####### get the eligible SNP:gene pair
cd /pi/manuel.garber-umw/human/skin/eQTLs/case_studies/DHX58/rs8081327
cat KRT_DESeq2_results_IFNvsPBS.txt | tr ' ' '\t' | awk '{if ($8<0.05) print $0}' | awk 'NR==1{print $1,$2,$3,$4,$5,$6,$7}NR>1{if ($4<-1.5 || $4>1.5) print $2,$3,$4,$5,$6,$7,$8}' > KRT_DESeq2_results_IFNvsPBS_sigDEgenes_padj0.05absLog2FC1.5.txt
cat KRT_DESeq2_results_IFNvsPBS_sigDEgenes_padj0.05absLog2FC1.5.txt | awk 'NR>1{print $1}' | sort | uniq > KRT_DEgenes.txt
cd ${dir}
cat ../KRT_DEgenes.txt | awk '{OFS=FS="\t"}{print "rs8081327",$0}' > SNP_gene_pairs.txt

####### 4. For every SNP-Gene pair, run the model with rank normalization and permutation test
# ^ note: for future runs, I should just run rank normalization in the first round of modeling to begin with. this round should just be the permutation test.
mkdir ${dir}/log/build_models
mkdir ${dir}/modelingResult
n_PCs=10
n_latentVar=10
rm ${dir}/commands.txt
while read pair;do
	snp=$(echo ${pair} | cut -d' ' -f1)
	gene=$(echo ${pair} | cut -d' ' -f2)
	echo "bash ${dir}/pipeline_DREG_KRT_perSNPGenePair_with_permutation.sh ${snp} ${gene} a${n_PCs}b${n_latentVar} ${dir} masteroutput_a${n_PCs}b${n_latentVar}" >> ${dir}/commands.txt
done < ${dir}/SNP_gene_pairs.txt
sleep 5
bash /pi/manuel.garber-umw/sshan/scripts/function_collapse_commands.sh ${dir} ${dir}/commands.txt ${dir}/commands.joined.txt 10
sleep 5
while read c; do
	echo ${c} | bsub -W 05:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=5000] -q long -e "${dir}/log/build_models_%J%I.err" -o "${dir}/log/build_models_%J%I.out"
done < ${dir}/commands.joined.txt
sleep 5
while [[ $(bjobs | grep ${jobname} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${jobname} | wc -l) "jobs remaining"; date;sleep 60; done


