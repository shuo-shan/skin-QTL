#!/bin/bash
# motif analysis for each promoter/enhancer cluster targeting each gene cluster
keyword=$1 # cluster1, cluster2, ..., cluster7
region=$2 # promoter or enhancer or pruned_promoter or filtered_enhancer
celltype=$3 # FRB, MEL, KRT

cd /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/HOMER/
mkdir ${keyword}_${region}_${celltype}
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/HOMER/${keyword}_${region}_${celltype}
cRE_bed=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_1000bp_flanking_window.bed

##### COMPILING BED FILES
# fetch gene cluster
#cd /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/heatmap_RNA/
cd /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/heatmap_RNA_round2/
this_gene_cluster=$(ls /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/heatmap_RNA_round2/ | grep reordered_ | grep .txt | grep log2FC1.5 | grep -v info | grep ${keyword}.txt)
cat ${this_gene_cluster} > ${dir}/temp.genes
cd ${dir}
# fetch corresponding gene enhancers
dict_gene_enhancer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_closest_gene_enhancer_links_allcts.txt
grep -f ${dir}/temp.genes ${dict_gene_enhancer} | grep -w ${celltype} | cut -f2 > temp.this_gene_cluster_enhancer
grep -f temp.this_gene_cluster_enhancer ${cRE_bed} | cut -f1-4 > this_gene_cluster_enhancer.bed
# clean-up
rm temp*

##### RUNNING HOMER
# TUTORIAL: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
outDir=${dir}
# INITIAL SET-UP:  generate configured genome file for hg38
#source activate fastQTL #--> change this to whatever conda environment you have.
#conda install -c bioconda homer # ver4.11
#perl /home/shuo.shan-umw/miniconda3/envs/fastQTL/share/homer/configureHomer.pl -install hg38 

# Run HOMER Motif enrichment and discovery for a given set of peaks
source activate fastQTL
PATH=$PATH:/home/shuo.shan-umw/miniconda3/envs/fastQTL/share/homer/bin/ #--> change this to conda install notes in the end.
pos_file=${dir}/this_gene_cluster_${region}.bed
findMotifsGenome.pl ${pos_file} hg38 ${outDir} -size given

