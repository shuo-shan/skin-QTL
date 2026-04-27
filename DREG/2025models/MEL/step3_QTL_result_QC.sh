#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 2:00
#BSUB -q short
#BSUB -J QTL.MEL[1-135]%70 # run 70 jobs concurrently
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/qc/qtl_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/qc/qtl_%J_%I.err"

sleep $((RANDOM % 10)) # reduces heavy I/O thundering-herd effect in parallel jobs
ct=MEL
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/results_QC

# Select chunk file safely
# Extract chunk base name and chunk index (000, 001, ...)
chunk=$(ls ${DIR}/chunks/pairs_chunk_*.tsv | sort -V | sed -n ${LSB_JOBINDEX}p)
chunk_base=$(basename "$chunk" .tsv)
chunk_id=$(echo "$chunk_base" | sed 's/pairs_chunk_//')

# -------- Subset genotype VCF file to 100K --------------------- #
echo "Performing QC for ${chunk}";date

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
bash -lc "
for cond in PBS IFNG IFNB TNF; do
  Rscript ${DIR}/step3_QTL_result_QC.R ${ct} \$cond eQTL  ${chunk_id}
done

for cond in IFNG IFNB TNF; do
  Rscript ${DIR}/step3_QTL_result_QC.R ${ct} \$cond reQTL ${chunk_id}
done
"
