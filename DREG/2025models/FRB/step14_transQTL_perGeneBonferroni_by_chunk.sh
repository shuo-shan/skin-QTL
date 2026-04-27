#!/bin/bash
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 1:00
#BSUB -q short
#BSUB -J QTL.FRB[1-101]%70 # run 70 jobs concurrently
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_bonferroni_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_bonferroni_%J_%I.err"

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs

ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/transQTL/bonferroni
mkdir -p ${DIR}/transQTL/bonferroni/PBS/eQTL
mkdir -p ${DIR}/transQTL/bonferroni/IFNG/eQTL ${DIR}/transQTL/bonferroni/IFNG/reQTL
mkdir -p ${DIR}/transQTL/bonferroni/IFNB/eQTL ${DIR}/transQTL/bonferroni/IFNB/reQTL
mkdir -p ${DIR}/transQTL/bonferroni/TNF/eQTL ${DIR}/transQTL/bonferroni/TNF/reQTL

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

## ---------------------------------------------------------- #
# Submit chunk to bonferroni Rscript, it calculates FWER adjusted p.val for gene:SNP pairs that passed QC
echo "Performing bonferroni for ${chunk_id}";date
cd ${DIR}

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step14_transQTL_perGeneBonferroni_by_chunk.R ${ct} ${cond} ${QTLtype} ${chunk_id}



