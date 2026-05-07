#!/bin/bash
#BSUB -J chrombpnet_score
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -q long
#BSUB -W 72:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step2_score.%J_%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step2_score.%J_%I.out"

# ============================================================
# DEFINE ONCE HERE — only touch these when changing runs
# ============================================================
BIAS_FACTOR=0.5
FOLD=0
sample=KRT_IFN
prefix=ATAC_${sample}_merged_alldonors

dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/${prefix}_biasFactor${BIAS_FACTOR}"
sif="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

# ============================================================
# SCRIPT TO SCORE VARIANTS
# ============================================================
rm -f ${dir}/run_inside_singularity_scoreSNPs.sh
cat > ${dir}/run_inside_singularity_scoreSNPs.sh << EOF
#!/bin/bash
module load deeptools
dir=${dir}
prefix=${prefix}
BIAS_FACTOR=${BIAS_FACTOR}
FOLD=${FOLD}

# run for each fold, then average with bigwigAverage (deeptools)
chrombpnet shap_scores \
    -m chrombpnet_nobias.h5 \
    -g hg38.fa \
    -r regions_of_interest.bed \
    -o shap_output/

echo "Pipeline complete"; date
EOF



echo "Waiting for filesystem sync..."; sleep 10
singularity exec --nv ${sif} bash ${dir}/run_inside_singularity.sh

