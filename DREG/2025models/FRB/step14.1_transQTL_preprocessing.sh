#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20000]"
#BSUB -W 72:00
#BSUB -q long
#BSUB -J transQTL.FRB
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_preprocessing_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_preprocessing_%J_%I.err"

ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/transQTL
mkdir -p ${DIR}/transQTL/eGene_QTL_pairs
mkdir -p ${DIR}/transQTL/QTL_tags
mkdir -p ${DIR}/transQTL/data

for cond in PBS IFNB IFNG TNF;do
        for QTLtype in eQTL reQTL;do
		if [[ "${cond}" == "PBS" && "${QTLtype}" == "reQTL" ]]; then
			continue
		fi
		mkdir -p ${DIR}/transQTL/chunks
		mkdir -p ${DIR}/transQTL/chunks/${cond}
		mkdir -p ${DIR}/transQTL/chunks/${cond}/${QTLtype}

		mkdir -p ${DIR}/transQTL/results
		mkdir -p ${DIR}/transQTL/results/${cond}
		mkdir -p ${DIR}/transQTL/results/${cond}/${QTLtype}
        done
done

# ---------------------------------- #
# make pair list between expressed/induced gene and lead cis-eQTLs
# ---------------------------------- #
# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

for cond in PBS IFNB IFNG TNF;do
	singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
		Rscript ${DIR}/step14.1_transQTL_preprocessing.R ${ct} ${cond} eQTL
done

for cond in IFNB IFNG TNF;do
	singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
		Rscript ${DIR}/step14.1_transQTL_preprocessing.R ${ct} ${cond} reQTL
done


# ---------------------------------- #
# compile genotype matrix for QTLs
# ---------------------------------- #
for cond in PBS IFNB IFNG TNF;do
	for QTLtype in eQTL reQTL;do
		if [[ "${cond}" == "PBS" && "${QTLtype}" == "reQTL" ]]; then
			continue
		fi
		bash ${DIR}/step14.1_transQTL_preprocessing_fetchGenotype.sh ${ct} ${cond} ${QTLtype}
	done
done
