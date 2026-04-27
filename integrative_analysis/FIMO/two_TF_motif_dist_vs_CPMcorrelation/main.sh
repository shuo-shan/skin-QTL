#!/bin/bash
# written by Crystal Shan 02/2024

TF1=GATA2
scriptDir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/two_TF_motif_dist_vs_CPMcorrelation

##### 1. For TF1, fetch the TF1-containing enhancers of TF1- highly-correlated and poorly-correlated genes.
cd /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/${TF1}
bash ${scriptDir}/function_fimo_fetch_region_with_motif_match_promoters_corrlated_vs_noncorrelated_part1.sh ${TF1}

##### 2. Get motif distance for each TF2 of interest.
cd /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/${TF1}
TFlist=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/heatmap_RNA_round2/TFs_in_TFcorrelationHeatmap.txt
cat ${TFlist} | grep -v ${TF1} > TF2_list.txt
while read TF2;do
	echo "bash ${scriptDir}/function_fimo_fetch_region_with_motif_match_promoters_corrlated_vs_noncorrelated_part2_perTF2.sh ${TF1} ${TF2}" >> commands.txt
done < TF2_list.txt
mkdir log
while read c; do
  echo ${c} | bsub -J ${TF1} -W 2:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -q long -e "./log/%J%I.err" -o "./log/%J%I.out"
done < commands.txt

while [[ $(bjobs | grep ${TF1} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${TF1} | wc -l) "jobs remaining"; sleep 5; done

###### clean-up
source activate fastQTL
pdfunite histograms_${TF1}_*.pdf merged_histograms.pdf
cat wilcoxon_test_pval_${TF1}_*.txt > wilcoxon_test_pval.txt
rm histograms_${TF1}_*.pdf wilcoxon_test_pval_${TF1}_*.txt
conda deactivate
