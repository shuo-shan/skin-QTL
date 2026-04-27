#!/bin/bash
# call Rscript to run susie on a gene of interest for each condition & QTLtype

# reduces heavy I/O thundering-herd effect in parallel jobs
sleep $((RANDOM % 20))

# set-up
ct=MEL
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
#mkdir -p ${DIR}/susie
#mkdir -p ${DIR}/susie/QC ${DIR}/susie/output_long
#mkdir -p ${DIR}/susie/PBS ${DIR}/susie/PBS/eQTL
#mkdir -p ${DIR}/susie/IFNG ${DIR}/susie/IFNG/eQTL ${DIR}/susie/IFNG/reQTL
#mkdir -p ${DIR}/susie/IFNB ${DIR}/susie/IFNB/eQTL ${DIR}/susie/IFNB/reQTL
#mkdir -p ${DIR}/susie/TNF ${DIR}/susie/TNF/eQTL ${DIR}/susie/TNF/reQTL

# ---------- fetch input arguments ------------------- #
g=$1   #g=ITGA1
echo "starting to find credible sets using SuSiE for ${g}"; date


# ------------ call Rscript to plot --------------- #
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step7_susie_perGene.R ${ct} ${g}
