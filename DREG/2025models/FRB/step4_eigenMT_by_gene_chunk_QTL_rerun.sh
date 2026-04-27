#!/bin/bash
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -W 24:00
#BSUB -q long
#BSUB -J chunk.FRB.rerun[18,35,105,183,184,230]%20
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/eigenMT_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/eigenMT_QTL_%J_%I.err"

# base memory is [mem=2000] for chunk[1-273]%100, after that some will fail. then also do chunk.rerun[18,35,105,183,184,230]%20 at mem=10000.

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs
ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/eigenMT
mkdir -p ${DIR}/eigenMT/PBS/eQTL
mkdir -p ${DIR}/eigenMT/IFNG/eQTL ${DIR}/eigenMT/IFNG/reQTL
mkdir -p ${DIR}/eigenMT/IFNB/eQTL ${DIR}/eigenMT/IFNB/reQTL
mkdir -p ${DIR}/eigenMT/TNF/eQTL ${DIR}/eigenMT/TNF/reQTL

# Convert LSB_JOBINDEX (1–120) → chunk_id (001–120)
chunk_id=$(printf "%03d" $((LSB_JOBINDEX)))

## ---------------------------------------------------------- #
# Submit chunk to eigenMT Rscript, it calculates FWER adjusted p.val for gene:SNP pairs that passed QC
echo "Performing eigenMT for ${chunk_id}";date
cd ${DIR}

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R FRB PBS eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R FRB IFNG eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R FRB IFNB eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R FRB TNF eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R FRB IFNG reQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R FRB IFNB reQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R FRB TNF reQTL ${chunk_id}







