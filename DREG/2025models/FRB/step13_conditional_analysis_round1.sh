#!/bin/bash
#BSUB -n 5
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 72:00
#BSUB -q long
#BSUB -J QTL.FRB[1-273]%70 # run 70 jobs concurrently
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/conditional1_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/conditional1_%J_%I.err"

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs
module load bcftools
module load htslib

ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/conditional_analysis_round1
mkdir -p ${DIR}/conditional_analysis_round1/results
mkdir -p ${DIR}/conditional_analysis_round1/results/PBS 
mkdir -p ${DIR}/conditional_analysis_round1/results/PBS/eQTL
mkdir -p ${DIR}/conditional_analysis_round1/results/IFNG
mkdir -p ${DIR}/conditional_analysis_round1/results/IFNG/eQTL
mkdir -p ${DIR}/conditional_analysis_round1/results/IFNG/reQTL
mkdir -p ${DIR}/conditional_analysis_round1/results/IFNB
mkdir -p ${DIR}/conditional_analysis_round1/results/IFNB/eQTL
mkdir -p ${DIR}/conditional_analysis_round1/results/IFNB/reQTL
mkdir -p ${DIR}/conditional_analysis_round1/results/TNF
mkdir -p ${DIR}/conditional_analysis_round1/results/TNF/eQTL
mkdir -p ${DIR}/conditional_analysis_round1/results/TNF/reQTL

# Select chunk file safely
chunk=$(ls ${DIR}/chunks/pairs_chunk_*.tsv | sort -V | sed -n ${LSB_JOBINDEX}p)
#chunk=${DIR}/chunks/pairs_chunk_096.tsv

# Extract chunk base name and chunk index (000, 001, ...)
chunk_base=$(basename "$chunk" .tsv)
chunk_id=$(echo "$chunk_base" | sed 's/pairs_chunk_//')

# -------- Set-up Singularity and Dependencies for Rscript ------ #
echo "Fitting QTL model for ${chunk}";date
cd ${DIR}

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step13_conditional_analysis_round1.R ${ct} 10 2 ${chunk_id} PBS eQTL
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step13_conditional_analysis_round1.R ${ct} 10 2 ${chunk_id} IFNB eQTL
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step13_conditional_analysis_round1.R ${ct} 10 2 ${chunk_id} IFNG eQTL
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step13_conditional_analysis_round1.R ${ct} 10 2 ${chunk_id} TNF eQTL
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step13_conditional_analysis_round1.R ${ct} 10 2 ${chunk_id} IFNB reQTL
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step13_conditional_analysis_round1.R ${ct} 10 2 ${chunk_id} IFNG reQTL
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step13_conditional_analysis_round1.R ${ct} 10 2 ${chunk_id} TNF reQTL
