#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=20400]
#BSUB -q long
#BSUB -W 08:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

# motif analysis for promoters or enhancers of a list of genes

# 1. change 'keyword' and 'genes' and 'region' for each run
keyword=melanocyte_expressed_genes_enhancer
genes=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/expressedGenes.txt
region=${dict_gene_enhancer}

# 2. set-up
cRE_bed=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_1000bp_flanking_window.bed
dict_gene_promoter=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_gene_promoter_links_allcts.txt
dict_gene_enhancer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_closest_gene_enhancer_links_allcts.txt
outDir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/HOMER/${keyword}
#cd /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/HOMER/
#mkdir ${outDir}; cd ${outDir}

# 3. create bed file
#bash /pi/manuel.garber-umw/sshan/scripts/function_rapid_fgrep.sh ${genes} ${region} temp1 rpdgrp 5000
#cat temp1 | grep -w MEL | cut -f2 > ${keyword}.txt
#bash /pi/manuel.garber-umw/sshan/scripts/function_rapid_fgrep.sh ${outDir}/${keyword}.txt ${cRE_bed} ${keyword}.bed rpdgrp 2000
pos_file=${outDir}/${keyword}.bed

# 4. Run HOMER Motif enrichment and discovery for a given set of peaks
# TUTORIAL: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
# INITIAL SET-UP:  generate configured genome file for hg38
#source activate fastQTL #--> change this to whatever conda environment you have.
#conda install -c bioconda homer # ver4.11
#perl /home/shuo.shan-umw/miniconda3/envs/fastQTL/share/homer/configureHomer.pl -install hg38
source activate fastQTL
PATH=$PATH:/home/shuo.shan-umw/miniconda3/envs/fastQTL/share/homer/bin/ #--> change this to conda install notes in the end.
findMotifsGenome.pl ${pos_file} hg38 ${outDir} -size given

rm -r rapid_fgrep_temp/
