#!/bin/bash
# run Rscript that performs coloc.susie or coloc.abf to compare PBS eQTL vs cytokine eQTL signal for any QTL gene for a cytokine across all 3 cytokines

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs

# ---------- fetch input arguments ------------------- #
celltype=$1
gene=$2
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${celltype}

# ---------- runs Rscript across all cytokines  ------------------- #
for cytokine in IFNG IFNB TNF;do
	export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
	singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step12_compare_PBSeQTL_cytokineeQTL_susie.coloc_perGene.R ${celltype} ${gene} ${cytokine}
	
	echo "wrote file to ${DIR}/coloc_susie/${cytokine}/coloc_susie_${gene}_coloc_summary.tsv"
	echo "wrote diagnostic plot to ${DIR}/coloc_susie/${cytokine}/coloc_susie_${gene}_coloc_diagnostic.pdf"
done
