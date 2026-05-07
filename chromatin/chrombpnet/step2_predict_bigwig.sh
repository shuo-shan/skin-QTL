#!/bin/bash
#BSUB -J chrombpnet_score
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -q gpu
#BSUB -W 72:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step2_score.%J_%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step2_score.%J_%I.out"

# ============================================================
# DEFINE ONCE HERE — only touch these when changing runs
# ============================================================
BIAS_FACTOR=0.5
sample=KRT_IFN
prefix=ATAC_${sample}_merged_alldonors

dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/${prefix}/biasFactor${BIAS_FACTOR}"
sif="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

module load deeptools

# ============================================================
# SCRIPT
# ============================================================
rm -f ${dir}/run_inside_singularity_postRun_analysis.sh
cat > ${dir}/run_inside_singularity_postRun_analysis.sh << EOF
#!/bin/bash
dir=${dir}
prefix=${prefix}
BIAS_FACTOR=${BIAS_FACTOR}

# run for each fold, then average with bigwigAverage (deeptools)
for fold_n in 0 1 2 3 4; do
	chrombpnet pred_bw \
	    -m ${dir}/fold${fold_n}/models/chrombpnet_nobias.h5 \
	    -g ${dir}/data/downloads/hg38.fa \
	    -c ${dir}/data/downloads/hg38.chrom.sizes \
	    -r ${dir}/data/downloads/peaks_no_blacklist.narrowPeak \ 
	    -o fold${fold_n}_pred/
done

# average across 5 folds
bigwigAverage \
    -b fold0.bw fold1.bw fold2.bw fold3.bw fold4.bw \
    -o averaged_chrombpnet.bw

# footprints (with average?)
chrombpnet footprints \
    -m chrombpnet_nobias.h5 \
    -g hg38.fa \
    -r peaks.bed \
    -motifs motifs.txt \
    -o footprints_output/


echo "Pipeline complete"; date
EOF



echo "Waiting for filesystem sync..."; sleep 10
singularity exec --nv ${sif} bash ${dir}/run_inside_singularity_postRun_analysis.sh


# output 
# For each variant, you get:
# log_counts_diff: log2(ALT/REF) predicted counts, Quantifies predicted accessibility change
# profile_jsd: JSD between REF and ALT predicted profiles, Detects profile shape changes
# Mean + std across models
