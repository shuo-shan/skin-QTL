#!/bin/bash
# motif analysis for each promoter/enhancer cluster targeting each gene cluster
cRE_file=$1 # full path, for example /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/heatmap_cRE/induced_enhancers_of_cts_induced_MEL.txt.

cd /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/HOMER/
keyword=$(echo ${cRE_file} | sed 's/.*\///g' | sed 's/\.txt//g')
mkdir ${keyword}
cd ${keyword}
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/HOMER/${keyword}
cRE_bed=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_300bp_flanking_window.bed

##### COMPILING BED FILES
grep -f ${cRE_file} ${cRE_bed} | cut -f1-4 > selected_enhancers.bed 

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
pos_file=${dir}/selected_enhancers.bed
findMotifsGenome.pl ${pos_file} hg38 ${outDir} -size given

