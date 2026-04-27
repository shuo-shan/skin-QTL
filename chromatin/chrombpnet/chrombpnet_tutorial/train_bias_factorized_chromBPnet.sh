#!/bin/bash
#BSUB -J chrombpnet
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10240]"
#BSUB -q gpu
#BSUB -W 720:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/trainBiasFactorizedChromBPnetModel.%J.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/trainBiasFactorizedChromBPnetModel.%J.out"


#### train ChromBPnet model using the bias model trained previously
#### activate singularity
singularity exec --nv /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif \
        bash -c '
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial
echo "training bias-factorized model using a pre-trained bias model now"; date
chrombpnet pipeline \
        -ibam ${dir}/data/downloads/merged.bam \
        -d "ATAC" \
        -g ${dir}/data/downloads/hg38.fa \
        -c ${dir}/data/downloads/hg38.chrom.sizes \
        -p ${dir}/data/downloads/peaks_no_blacklist.bed \
        -n ${dir}/data/output_negatives.bed \
        -fl ${dir}/data/splits/fold_0.json \
	-b ${dir}/bias_model/models/models/k562_bias.h5 \
        -o ${dir}/chrombpnet_model/

'
