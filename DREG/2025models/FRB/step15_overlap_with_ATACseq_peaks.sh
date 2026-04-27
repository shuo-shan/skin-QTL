#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 08:00
#BSUB -J overlapATAC.FRB
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/overlapATAC_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/overlapATAC_%J.err"

# overlap sig eQTLs/reQTLs with ATACseq peaks
# version1: phase1 ATACseq file F25,F49,F55, mapped with in-house pipeline, PBS and IFNG samples only, 1kb flanking region around summit
# version2: phase1+2 ATACseq file F25,F49,F55,F108,F109,F110, mapped with encode pipeline, PBS and IFNG samples only, 1kb flanking region around summit

# set-up
ct=FRB
cond=IFNG
QTLtype=eQTL

DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/overlapATAC
cd ${DIR}

# ------- compile significant eQTL table across all eGenes of this condition x QTLtype ---------
file2="${DIR}/eigenMT/results/${ct}_${cond}_${QTLtype}.eigenMT.txt"
outF="${DIR}/overlapATAC/QTL_stats_${ct}_${cond}_${QTLtype}.txt"
rm -f "${outF}"

first=1
for file1 in "${DIR}"/results_QC/modeling_stats_postQC_${ct}_${cond}_${QTLtype}_*.txt; do
    echo "${file1}"

    awk '
    BEGIN { FS=OFS="\t" }

    NR==FNR {
        if (FNR==1) {
            for (i=1; i<=NF; i++) {
                if ($i=="gene") gene_col=i
                if ($i=="Meff") meff_col=i
                if ($i=="q_gene") q_col=i
            }
            next
        }
        if (($(q_col)+0) < 0.05) {
            meff[$gene_col] = $(meff_col)+0
            cutoff[$gene_col] = 0.05 / ($(meff_col)+0)
        }
        next
    }

    FNR==1 {
        for (i=1; i<=NF; i++) {
            if ($i=="gene") gene2_col=i
            if ($i=="p") p_col=i
        }
        if (print_header) print $0, "Meff", "p_nominal_cutoff"
        next
    }

    ($gene2_col in meff) && (($(p_col)+0) <= cutoff[$gene2_col]) {
        print $0, meff[$gene2_col], cutoff[$gene2_col]
    }
    ' print_header=$first "${file2}" "${file1}" >> "${outF}"

    first=0
done

# ------- annotate the CHR and POS for QTLs -----------

# ------- join with ATACseq peaks --------- #
merged_peaks=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_peaks.narrowPeak
dynamics_FRB=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_FRB.txt
dynamics_MEL=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_MEL.txt
dynamics_KRT=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_KRT.txt


