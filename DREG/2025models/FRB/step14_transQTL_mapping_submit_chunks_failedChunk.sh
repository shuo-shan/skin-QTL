#!/bin/bash
#BSUB -n 5
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 72:00
#BSUB -q long
#BSUB -J QTL.FRB
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_%J_%I.err"

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs
module load bcftools
module load htslib

ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}

cond=PBS
QTLtype=eQTL
chunk_id=028
chunk=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/transQTL/chunks/PBS/eQTL/eGene_QTL_pairs_chunk_028.tsv

echo "Running task:"
echo "  ct=${ct}"
echo "  cond=${cond}"
echo "  QTLtype=${QTLtype}"
echo "  chunk_id=${chunk_id}"
echo "  chunk=${chunk}"
date

cd "${DIR}"

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
    Rscript "${DIR}/step14_transQTL_run_QTL_chunk.R" "${ct}" 10 2 "${chunk_id}" "${cond}" "${QTLtype}"
