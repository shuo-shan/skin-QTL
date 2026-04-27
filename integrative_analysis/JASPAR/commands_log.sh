# log of jaspar enrichment analysis commands

# 04/28/2023
jaspar_dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/JASPAR/jaspar_enrichment/bin
keyword=cluster10_enhancers
genes=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/heatmap_RNA_round2/reordered_inducedGenes_log2FC1.5_cluster10.txt
region=${dict_gene_enhancer}

cRE_bed=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_1000bp_flanking_window.bed
dict_gene_promoter=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_gene_promoter_links_allcts.txt
dict_gene_enhancer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_closest_gene_enhancer_links_allcts.txt
outDir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/JASPAR/${keyword}



${jaspar_dir}/JASPAR_enrich.sh oneSetBg ${jaspar_dir} 


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

jaspar_dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/JASPAR/jaspar_enrichment/bin
${jaspar_dir}/JASPAR_enrich.sh oneSetBg <loladb_dir> <S bed> <U bed> <output dir> <API_URL> <n_cores>
# <loladb_dir> is the directory containing all the LOLA database batches for an assembly
# <API_URL> is the URL to the matrix API in JASPAR, e.g. http://jaspar.genereg.net/api/v1/matrix/
# <n_cores> is the number of cores to use when parallelizing.
# <S bed> is the genomic regions we are interested in, .bed format
# <U bed> is the universie set of genomic regions. every region in S should overlap with one region in U.

# example 1:
# this gives us a list of enriched TFs in the regions. You can make a beeswarm plot with the -log10(P-value) and color code the TF class
${jaspar_dir}/JASPAR_enrich.sh oneSetBg ${hg38_dir} promoters.bed all_regions.bed ${output_dir}



# Differential Enrichment
${jaspar_dir}/JASPAR_enrich.sh twoSets <loladb_dir> <S1 bed> <S2 bed> <output dir> <API URL> <n_cores>
# the <U bed> is not needed here because JASPAR treats the U as the merge between S1 and S2 bed files.

