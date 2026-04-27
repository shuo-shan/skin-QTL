#!/bin/bash
#BSUB -n 5
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 24:00
#BSUB -q long
#BSUB -J QTL.FRB.mapping[105]%160
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_%J_%I.err"

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs
module load bcftools
module load htslib

ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
taskfile=${DIR}/transQTL/task_manifest.tsv

# read one line from manifest
task=$(sed -n "${LSB_JOBINDEX}p" "${taskfile}")

if [[ -z "${task}" ]]; then
    echo "No task found for LSB_JOBINDEX=${LSB_JOBINDEX}"
    exit 1
fi

cond=$(echo "${task}" | cut -f1)
QTLtype=$(echo "${task}" | cut -f2)
chunk_id=$(echo "${task}" | cut -f3)
chunk=$(echo "${task}" | cut -f4)

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
    Rscript "${DIR}/step14.3_transQTL_run_QTL_chunk.R" "${ct}" 10 2 "${chunk_id}" "${cond}" "${QTLtype}"
