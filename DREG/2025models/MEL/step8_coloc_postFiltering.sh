#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=1000]
#BSUB -q long
#BSUB -W 08:00
#BSUB -J job_submitter
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/job_submitter_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/job_submitter_%J.err"
# batch submit jobs to run SuSiE on QTL genes (fdr 0.05)

# set-up
ct=MEL
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}

cd ${DIR}

# after jobs run, summarize and sort by PP.H4 descending.
QTLtype=eQTL
for cond in PBS IFNG IFNB TNF;do
	echo "processing ${cond} ${QTLtype}"
	cd ${DIR}/coloc/vitiligo/${cond}/${QTLtype}

	outF=${DIR}/coloc/vitiligo/output_long/coloc_vitiligo_${cond}_${QTLtype}
	head -1 coloc_summary_A1BG-AS1.tsv > ${outF} 
	for f in coloc_summary_*.tsv; do
		awk 'NR>1' ${f}  >> ${outF}
	done

	outFsorted=${outF}.sorted.txt
	{
		head -1 ${outF}
		tail -n +2 ${outF} | sort -t$'\t' -k10,10gr
	} > ${outFsorted}
	rm ${outF}
done

QTLtype=reQTL
for cond in IFNG IFNB TNF;do
	echo "processing ${cond} ${QTLtype}"
	cd ${DIR}/coloc/vitiligo/${cond}/${QTLtype}

	outF=${DIR}/coloc/vitiligo/output_long/coloc_vitiligo_${cond}_${QTLtype}
	head -1 coloc_summary_A1BG-AS1.tsv > ${outF} 
	for f in coloc_summary_*.tsv; do
		awk 'NR>1' ${f}  >> ${outF}
	done

	outFsorted=${outF}.sorted.txt
	{
		head -1 ${outF}
		tail -n +2 ${outF} | sort -t$'\t' -k10,10gr
	} > ${outFsorted}
	rm ${outF}
done

# Pick SNP:gene pairs where PPH4 > 70%
outF=${DIR}/coloc/vitiligo/output_long/coloc_vitiligo_PPH4_70.txt
head -1 ${DIR}/coloc/vitiligo/output_long/coloc_vitiligo_IFNG_eQTL.sorted.txt > ${outF}
for QTLtype in eQTL reQTL; do
	for cond in PBS IFNG IFNB TNF; do
		inFbase=${DIR}/coloc/vitiligo/output_long/coloc_vitiligo_${cond}_${QTLtype}
		inFsorted=${inFbase}.sorted.txt

		if [[ -f ${inFsorted} ]]; then
			echo "processing ${inFsorted}"
			awk 'NR>1{if ($10 > 0.70) print $0}' ${inFsorted} >> ${outF}
		fi
	done
done

# move top gene plots to folder
mkdir -p ${DIR}/coloc/vitiligo/output_long/PPH4_70
awk 'NR>1{print $2}' ${outF} | sort -u > ${DIR}/coloc/vitiligo/output_long/top_genes.txt
while read g;do
	echo ${g}
	cp ${DIR}/coloc/vitiligo/plots/MEL_${g}.locus_tracks.pdf ${DIR}/coloc/vitiligo/output_long/PPH4_70
done < ${DIR}/coloc/vitiligo/output_long/top_genes.txt


#### next I want to know distinct causal variants PPH3 > 70% ---- ####
# Pick SNP:gene pairs where PPH4 > 70%
outF=${DIR}/coloc/vitiligo/output_long/coloc_vitiligo_PPH3_70.txt
head -1 ${DIR}/coloc/vitiligo/output_long/coloc_vitiligo_IFNG_eQTL.sorted.txt > ${outF}
for QTLtype in eQTL reQTL; do
        for cond in PBS IFNG IFNB TNF; do
                inFbase=${DIR}/coloc/vitiligo/output_long/coloc_vitiligo_${cond}_${QTLtype}
                inFsorted=${inFbase}.sorted.txt

                if [[ -f ${inFsorted} ]]; then
                        echo "processing ${inFsorted}"
                        awk 'NR>1{if ($9 > 0.70) print $0}' ${inFsorted} >> ${outF}
                fi
        done
done

mkdir -p ${DIR}/coloc/vitiligo/output_long/PPH3_70
awk 'NR>1{print $2}' ${outF} | sort -u > ${DIR}/coloc/vitiligo/output_long/top_genes.txt
while read g;do
        echo ${g}
        cp ${DIR}/coloc/vitiligo/plots/MEL_${g}.locus_tracks.pdf ${DIR}/coloc/vitiligo/output_long/PPH3_70
done < ${DIR}/coloc/vitiligo/output_long/top_genes.txt
rm ${DIR}/coloc/vitiligo/output_long/top_genes.txt
