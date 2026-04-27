#!/bin/bash
#BSUB -J trainbias
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10240]"
#BSUB -q gpu
#BSUB -W 720:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/MEL_IFNG/log/trainBiasModel.%J.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/MEL_IFNG/log/trainBiasModel.%J.out"


#### activate singularity
singularity exec --nv /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif \
	bash -c '
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/MEL_IFNG
mkdir -p ${dir}/bias_model
mkdir -p ${dir}/data/output
prefix=ATAC_F108_MEL_IFNG_S1_100Ksample

### this is how to train a bias model
echo "training bias model now"; date
chrombpnet bias pipeline \
        -ibam ${dir}/data/downloads/merged.bam \
        -d "ATAC" \
        -g ${dir}/data/downloads/hg38.fa \
        -c ${dir}/data/downloads/hg38.chrom.sizes \
        -p ${dir}/data/downloads/peaks_no_blacklist.bed \
        -n ${dir}/data/output_negatives.bed \
        -fl ${dir}/data/splits/fold_0.json \
        -b 0.5 \
        -o ${dir}/bias_model/ \
	-s 42 \
        -fp ${prefix}_fold0 
echo "done";date
'

