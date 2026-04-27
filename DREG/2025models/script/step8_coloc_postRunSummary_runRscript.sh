#!/bin/bash
# call Rscript to run coloc on a gene of interest for each condition & QTLtype

# reduces heavy I/O thundering-herd effect in parallel jobs
sleep $((RANDOM % 20))

# ---------- fetch input arguments ------------------- #
ct=$1
g=$2   #g=ITGA1
trait=$3
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
echo "starting to perform coloc on all GWAS traits for ${g}"; date

# ------------ call Rscript to plot --------------- #
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
echo ${trait}; date
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif Rscript ${DIR}/step8_coloc_join_coloc_summary_by_stats_perGene.R ${ct} ${g} ${trait}
