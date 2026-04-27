#!/bin/bash
# perform BH correction on locally adjusted pvalues across all genes
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -W 2:00
#BSUB -q long
#BSUB -J BH
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/eigenMT_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/eigenMT_QTL_%J_%I.err"

ct=MEL
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/eigenMT/results

# concatenate all eigenMT result tables 
for cond in PBS IFNG IFNB TNF;do
	echo "concatenating eigenMT results for ${cond} eQTLs"
	cd ${DIR}/eigenMT/${cond}/eQTL
	head -1 ${ct}_${cond}_eQTL_000.eigenMT.tsv > ${DIR}/eigenMT/results/${ct}_${cond}_eQTL.eigenMT.txt
	for f in ${ct}_${cond}_eQTL*.tsv;do
		awk 'NR>1' ${f} >> ${DIR}/eigenMT/results/${ct}_${cond}_eQTL.eigenMT.txt
	done
done

for cond in IFNG IFNB TNF;do
        echo "doing BH correction for ${cond} reQTL across genes"
	echo "concatenating eigenMT results for ${cond} reQTLs"
        cd ${DIR}/eigenMT/${cond}/reQTL
        head -1 ${ct}_${cond}_reQTL_000.eigenMT.tsv > ${DIR}/eigenMT/results/${ct}_${cond}_reQTL.eigenMT.txt
        for f in ${ct}_${cond}_reQTL*.tsv;do
                awk 'NR>1' ${f} >> ${DIR}/eigenMT/results/${ct}_${cond}_reQTL.eigenMT.txt
        done
done


# perform BH
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

for cond in PBS IFNG IFNB TNF;do
	echo "doing BH correction for ${cond} eQTL across genes"
	cd ${DIR}/eigenMT/results
	singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
        Rscript ${DIR}/step5_BHcorrection.R ${ct}_${cond}_eQTL.eigenMT.txt
done

for cond in IFNG IFNB TNF;do
	echo "doing BH correction for ${cond} reQTL across genes"
	cd ${DIR}/eigenMT/results
	singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
        Rscript ${DIR}/step5_BHcorrection.R ${ct}_${cond}_reQTL.eigenMT.txt
done
