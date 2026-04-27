#!/bin/bash
#BSUB -J predict_bw
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=9240]"
#BSUB -q gpu
#BSUB -W 720:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/KRT_PBS_dedupped_bias_threshold_factor_0.7/log/predictbw.%J.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/KRT_PBS_dedupped_bias_threshold_factor_0.7/log/predictbw.%J.out"

## interactive gpu session
# start the gpu session
bsub -Is -q gpu -W 8:00 -R rusage[mem=100G] /bin/bash
# start the singularity image
singularity shell --nv /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif

dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/KRT_PBS_dedupped_bias_threshold_factor_0.7"
cd ${dir}
mkdir -p bigwig
mkdir -p footprints
prefix=ATAC_F111_KRT_PBS_S1_fold0
BIAS_MODEL=${dir}/bias_model/models/${prefix}_bias.h5
CHROMBPNET_MODEL=${dir}/chrombpnet_model/models/chrombpnet.h5
CHROMBPNET_MODEL_NB=${dir}/chrombpnet_model/models/chrombpnet_nobias.h5
GENOME=${dir}/data/downloads/hg38.fa
CHROM_SIZES=${dir}/data/downloads/hg38.chrom.sizes
OUTPUT_PREFIX=${prefix}

# make 10-column bed file, need at least 2 rows
# Make sure that regions in the input bed file can be expanded to inputlen (default to 2114) regions without overflowing out of the chromosomes.
# DHX58: chr17:42,106,878-42,113,433
# IRF4 intron: chr6:402,937-404,100
#echo -e "chr6\t402443\t404557\tIRF4_intron\t0\t.\t0\t0\t0\t-1" > ${dir}/data/downloads/region.bed
#echo -e "chr17\t42104878\t42106992\tDHX58\t0\t.\t0\t0\t0\t-1" >> ${dir}/data/downloads/region.bed
echo -e "chr17\t42107360\t42109474\tregion_1\t0\t.\t0\t0\t0\t-1
chr17\t42109474\t42111588\tregion_2\t0\t.\t0\t0\t0\t-1
chr17\t42111588\t42113702\tregion_3\t0\t.\t0\t0\t0\t-1
chr17\t42113702\t42115816\tregion_4\t0\t.\t0\t0\t0\t-1
chr17\t42115816\t42117930\tregion_5\t0\t.\t0\t0\t0\t-1
chr17\t42117930\t42120044\tregion_6\t0\t.\t0\t0\t0\t-1
chr17\t42120044\t42122158\tregion_7\t0\t.\t0\t0\t0\t-1
chr17\t42122158\t42124272\tregion_8\t0\t.\t0\t0\t0\t-1
chr17\t42124272\t42126386\tregion_9\t0\t.\t0\t0\t0\t-1
chr17\t42126386\t42128500\tregion_10\t0\t.\t0\t0\t0\t-1" > ${dir}/data/downloads/region.bed
REGIONS=${dir}/data/downloads/region.bed

### ---------------- run pred_bw
chrombpnet pred_bw -bm ${BIAS_MODEL} -cm ${CHROMBPNET_MODEL} -cmb ${CHROMBPNET_MODEL_NB} -r ${REGIONS} -g ${GENOME} -c ${CHROM_SIZES} -op ${OUTPUT_PREFIX}

### ---------------- generate contribution score bigwigs
chrombpnet contribs_bw -m ${CHROMBPNET_MODEL_NB} -r ${REGIONS} -g ${GENOME} -c ${CHROM_SIZES} -op ${OUTPUT_PREFIX}

### ---------------- TF marginal footprinting
CHR_FOLD_PATH=${dir}/data/splits/fold_0.json
MOTIFS_TO_PWM=${dir}/data/downloads/pwm.tsv 
OUTPUT_PREFIX=${prefix}.marginalfp
# generate non-peak regions
shuf -n 100000 --random-source=<(yes 42) ${dir}/data/output_negatives.bed > ${dir}/data/downloads/nonpeaks_100ksample.bed
NONPEAK_REGIONS=${dir}/data/downloads/nonpeaks_100ksample.bed
# The argument -pwm_f is a path to a TSV file containing motifs in first column (e.g. Tn5) and motif string (e.g. GCACAGTACAGAGCTG) to use for footprinting in second column. A default file is provided in the data folder for reference (https://github.com/kundajelab/chrombpnet/blob/master/chrombpnet/data/motif_to_pwm.TF.tsv)
chrombpnet footprints -m ${CHROMBPNET_MODEL} -r ${NONPEAK_REGIONS} -g ${GENOME} -fl ${CHR_FOLD_PATH} -op ${OUTPUT_PREFIX} -pwm_f ${MOTIFS_TO_PWM} 


### ---------------- variant effect prediction
SNP_DATA=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/temp_snps.tsv
OUTPUT_PREFIX=${prefix}.snpeffect
chrombpnet snp_score -snps ${SNP_DATA} -m ${CHROMBPNET_MODEL} -g ${GENOME} -op ${OUTPUT_PREFIX}
