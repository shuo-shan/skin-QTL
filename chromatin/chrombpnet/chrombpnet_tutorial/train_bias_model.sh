#!/bin/bash
#BSUB -J chrombpnet
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10240]"
#BSUB -q gpu
#BSUB -W 720:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/trainBiasModel.%J.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/trainBiasModel.%J.out"

#### start by downloading a pre-trained bias model
#mkdir -p ${dir}/bias_model
#wget https://storage.googleapis.com/chrombpnet_data/input_files/bias_models/ATAC/ENCSR868FGK_bias_fold_0.h5 -O ${dir}/bias_model/ENCSR868FGK_bias_fold_0.h5

### train a bias-factorized ChromBPNet model using a pre-trained bias model. 
# I'll learn how to train a bias model next
#### activate singularity
singularity exec --nv /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif \
	bash -c '
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial

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
        -o ${dir}/bias_model/models \
        -fp k562 
echo "done";date
'

