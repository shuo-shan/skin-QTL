#!/bin/bash
# perform BH correction on locally adjusted pvalues across all genes
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -W 02:00
#BSUB -q long
#BSUB -J BH
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/BH_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/BH_QTL_%J_%I.err"

ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
OUTDIR=${DIR}/transQTL/resultsBHcorrected
mkdir -p ${OUTDIR}

# concatenate all results
for cond in PBS IFNG IFNB TNF;do
	for QTLtype in eQTL reQTL;do
		if [[ "$cond" == "PBS" && "$QTLtype" == "reQTL" ]];then
			continue
		fi
		echo "concatenating trans-QTL mapping results for ${cond} ${QTLtype}"
		cd ${DIR}/transQTL/results/${cond}/${QTLtype}
		head -1 result_001.tsv > ${OUTDIR}/${ct}_${cond}_${QTLtype}.result.txt
		for f in *.tsv;do
			awk 'NR>1' ${f} >> ${OUTDIR}/${ct}_${cond}_${QTLtype}.result.txt
		done
	done
done

# perform BH across genes using the pmin of bonferroni-corrected pvalue
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

for cond in PBS IFNG IFNB TNF;do
	for QTLtype in eQTL reQTL;do
                if [[ "$cond" == "PBS" && "$QTLtype" == "reQTL" ]];then
                        continue
                fi
		echo "doing BH correction for ${cond} ${QTLtype} across all pairs"
		cd ${OUTDIR}
		singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
		Rscript ${DIR}/step14.4_transQTL_flatBH.R ${ct}_${cond}_${QTLtype}.result.txt
	done
done

