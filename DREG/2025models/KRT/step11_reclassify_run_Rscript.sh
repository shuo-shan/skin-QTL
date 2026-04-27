#!/bin/bash
# run Rscript that reclassifies any QTL gene for a cytokine across all 3 cytokines

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs

# ---------- fetch input arguments ------------------- #
celltype=$1
gene=$2
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${celltype}

# ---------- runs Rscript across all cytokines  ------------------- #
for cytokine in IFNG IFNB TNF;do
	export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
	singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step10_reclassify_QTL_genes_${celltype}.R ${celltype} ${cytokine} ${gene}
	
	echo "wrote file to ${DIR}/reclassified/${cytokine}/reclassified_${gene}.txt"
done
