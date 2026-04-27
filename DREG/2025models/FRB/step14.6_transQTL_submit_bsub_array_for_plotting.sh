#!/bin/bash
#BSUB -J transQTL[1-16]%270 #change this number based on number of SNPs
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -W 1:00
#BSUB -q short
#BSUB -R "rusage[mem=1500]"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_plot_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_plot_%J_%I.err"

ct=FRB
condition=TNF
QTLtype=eQTL

# SNP list (one per line, matching array indices 1-21)
snps_file=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}/transQTL/snps.txt

# Read SNP for this array index (sed is 1-indexed, matching LSB_JOBINDEX)
SNP=$(sed -n "${LSB_JOBINDEX}p" ${snps_file})


echo "working on ${SNP}";date
bash step14.6_transQTL_compile_result_for_snp.sh ${SNP} ${ct} ${condition} ${QTLtype}
