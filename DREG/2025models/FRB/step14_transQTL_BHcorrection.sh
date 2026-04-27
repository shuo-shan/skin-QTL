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
mkdir -p ${DIR}/transQTL/bonferroni/results

# concatenate all bonferroni result tables 
for cond in PBS IFNB IFNG TNF;do
        for QTLtype in eQTL reQTL;do
                if [[ "${cond}" == "PBS" && "${QTLtype}" == "reQTL" ]]; then
                        continue
                fi
		echo "concatenating bonferroni results for ${cond} ${QTLtype}"
		cd ${DIR}/transQTL/bonferroni/${cond}/${QTLtype}
		head -1 ${ct}_${cond}_${QTLtype}_001.bonferroni.tsv > ${DIR}/transQTL/bonferroni/results/${ct}_${cond}_${QTLtype}.bonferroni.txt
		for f in ${ct}_${cond}_${QTLtype}*.tsv;do
			awk 'NR>1' ${f} >> ${DIR}/transQTL/bonferroni/results/${ct}_${cond}_${QTLtype}.bonferroni.txt
		done
	done
done


# perform BH
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

for cond in PBS IFNB IFNG TNF;do
        for QTLtype in eQTL reQTL;do
                if [[ "${cond}" == "PBS" && "${QTLtype}" == "reQTL" ]]; then
                        continue
                fi
		echo "doing BH correction for ${cond} ${QTLtype} across genes"
		cd ${DIR}/transQTL/bonferroni/results
		singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
		Rscript ${DIR}/step14_transQTL_BHcorrection.R ${ct}_${cond}_${QTLtype}.bonferroni.txt
	done
done

