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
celltype=FRB
jobname=${celltype}
dir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/reQTL_betacomparison/FRB
cd ${dir}
mkdir plot log

# copy the covariate files, expressed gene files from basic DREG pipeline
olddir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/FRB_new
cp ${olddir}/expressed_genes_tss* ${dir}
cp ${olddir}/snps_near_expressed_genes.bed ${dir}

# get the eligible SNP:gene pair (405,116)
cat ${olddir}/masteroutput_all_with_colnames.txt | awk '{OFS="\t"}NR>1{if ( ($7!="." && $7<0.01) || ($15!="." && $15<0.01) || ($23!="." && $23<0.01) ) print $1"_"$2}' | sort | uniq | tr '_' ' ' > snp_gene_pairs.txt

# for each SNP:gene pair, run the beta-comparison model with 10,000 permutations
mkdir ${dir}/log/build_models
mkdir ${dir}/modelingResult
cat ${dir}/snp_gene_pairs.txt | awk -v dir=${dir} '{print "bash "dir"/pipeline_betacomparison_FRB_perSNPGenePair.sh "$0" "dir" masteroutput"}' > ${dir}/commands.txt
sleep 3
bash /pi/manuel.garber-umw/sshan/scripts/function_collapse_commands.sh ${dir} ${dir}/commands.txt ${dir}/commands.joined.txt 1000
sleep 5
while read c; do
        echo ${c} | bsub -W 72:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=800] -q long -e "${dir}/log/build_models_%J%I.err" -o "${dir}/log/build_models_%J%I.out"
done < ${dir}/commands.joined.txt

while [[ $(bjobs | grep ${jobname} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${jobname} | wc -l) "jobs remaining"; date;sleep 60; done



