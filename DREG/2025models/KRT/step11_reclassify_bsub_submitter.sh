#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=1000]
#BSUB -q long
#BSUB -W 08:00
#BSUB -J job_submitter
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/KRT/logs/job_submitter_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/KRT/logs/job_submitter_%J.err"
# batch submit jobs to run SuSiE on QTL genes (fdr 0.05)

# set-up
ct=KRT
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/reclassified
LOG_ROOT=${DIR}/logs/reclassify; mkdir -p ${LOG_ROOT}

gene_file=${DIR}/data/fdr05_genelist.txt
cat ${DIR}/eigenMT/results/*fdr05_genelist.txt | sort -u | sed '/^$/d' > ${gene_file}

# compile commands
CMD=${DIR}/step10_reclassify_commands.txt
rm -f ${CMD}
while read gene;do
	echo "bash ${DIR}/step10_reclassify_run_Rscript.sh ${ct} ${gene}" >> ${CMD}
done < ${gene_file}

# create helper file
cat > ${DIR}/run_task.sh <<'EOF'
#!/bin/bash
set -euo pipefail
cmd=$(sed -n "${LSB_JOBINDEX}p" "$1"); echo "RUN: $cmd"; eval "$cmd"
EOF
chmod +x ${DIR}/run_task.sh

# Submit BSUB jobs
MAX_CONCURRENT=200
N=$(wc -l < ${CMD})
bsub -q short -W 0:30 -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" \
  -J "reclassify.${ct}[1-${N}]%${MAX_CONCURRENT}" \
  -o "${LOG_ROOT}/reclassify_%J_%I.out" \
  -e "${LOG_ROOT}/reclassify_%J_%I.err" \
  ${DIR}/run_task.sh ${CMD}


# Monitor jobs' progress
while [[ $(bjobs -w | grep reclassify.${ct} | wc -l) != 0 ]] ; do echo $(bjobs -w | grep reclassify.${ct} | wc -l) "jobs remaining"; sleep 5; done

# Check out if any jobs failed
cd ${LOG_ROOT}
rm -f ${LOG_ROOT}/temp.failedjobs.txt
for f in *out; do
	if ! grep -q "Successfully completed" ${f}; then
		echo "${f}" >> "${LOG_ROOT}/temp.failedjobs.txt"
	fi
done
echo "Total jobs: $(ls *.out | wc -l)  |    Failed jobs: $(wc -l temp.failedjobs.txt)"

# After jobs run and QC'ed, summarize
for cytokine in IFNB IFNG TNF;do
	summary=${DIR}/reclassified/reclassify_summary_${cytokine}.txt
	rm -f ${summary}

	cd ${DIR}/reclassified/${cytokine}
	this_file=$(ls ${DIR}/reclassified/${cytokine} | head -1)
	awk '{OFS=FS="\t"}NR==1{print $0}' ${this_file} > ${summary}
	for f in reclassified_*.txt; do
		awk '{OFS=FS="\t"}NR>1{print $0}' ${f} >> ${summary}
	done
	echo "${summary}"
done

Rscript ${DIR}/step10_reclassify_postRunCollection.R






