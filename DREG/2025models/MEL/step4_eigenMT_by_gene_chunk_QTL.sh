#!/bin/bash
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=5000]"
#BSUB -W 72:00
#BSUB -q long
#BSUB -J chunk.MEL[1-120]%120 # run 120 jobs concurrently
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/eigenMT_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/eigenMT_QTL_%J_%I.err"

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs
module load bcftools
ct=MEL
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/eigenMT
mkdir -p ${DIR}/eigenMT/PBS/eQTL
mkdir -p ${DIR}/eigenMT/IFNG/eQTL ${DIR}/eigenMT/IFNG/reQTL
mkdir -p ${DIR}/eigenMT/IFNB/eQTL ${DIR}/eigenMT/IFNB/reQTL
mkdir -p ${DIR}/eigenMT/TNF/eQTL ${DIR}/eigenMT/TNF/reQTL

all_genes=${DIR}/results_QC/all_master_pairs_${ct}_genes.txt
all_pairs=${DIR}/results_QC/all_master_pairs_${ct}.txt
all_modelstats=${DIR}/permutation/model_stats/all_model_results.txt
vcf="/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz"

# Convert LSB_JOBINDEX (1–120) → chunk_id (000–119)
chunk_id=$(printf "%03d" $((LSB_JOBINDEX - 1)))

## ---------------------------------------------------------- #
# Submit chunk to eigenMT Rscript
echo "Performing eigenMT for ${chunk_id}";date
cd ${DIR}

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R MEL PBS eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R MEL IFNG eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R MEL IFNB eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R MEL TNF eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R MEL IFNG reQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R MEL IFNB reQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_eigenMT_gene_QTL_perGeneChunk.R MEL TNF reQTL ${chunk_id}







