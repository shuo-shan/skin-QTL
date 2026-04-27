#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=1000]
#BSUB -q long
#BSUB -W 08:00
#BSUB -J job_submitter
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/job_submitter_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/job_submitter_%J.err"
# batch submit jobs to run SuSiE on QTL genes (fdr 0.05)

# set-up
ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/coloc_susie
LOG_ROOT=${DIR}/logs/coloc_susie; mkdir -p ${LOG_ROOT}

gene_file=${DIR}/data/fdr05_genelist.txt
cat ${DIR}/eigenMT/results/*fdr05_genelist.txt | sort -u | sed '/^$/d' > ${gene_file}

# compile commands
CMD=${DIR}/step12_commands.txt
rm -f ${CMD}
while read gene;do
	echo "bash ${DIR}/step12_compare_PBSeQTL_cytokineeQTL_susie.coloc_run_Rscript.sh ${ct} ${gene}" >> ${CMD}
done < ${gene_file}

# create helper file
cat > ${DIR}/run_task.sh <<'EOF'
#!/bin/bash
set -euo pipefail
cmd=$(sed -n "${LSB_JOBINDEX}p" "$1"); echo "RUN: $cmd"; eval "$cmd"
EOF
chmod +x ${DIR}/run_task.sh

# Submit BSUB jobs
MAX_CONCURRENT=150
N=$(wc -l < ${CMD})
bsub -q short -W 0:59 -n 1 -R "span[hosts=1]" -R "rusage[mem=3000]" \
  -J "PBSvsCyto.${ct}[1-${N}]%${MAX_CONCURRENT}" \
  -o "${LOG_ROOT}/PBSvsCyto_%J_%I.out" \
  -e "${LOG_ROOT}/PBSvsCyto_%J_%I.err" \
  ${DIR}/run_task.sh ${CMD}


# Monitor jobs' progress
while [[ $(bjobs -w | grep PBSvsCyto.${ct} | wc -l) != 0 ]] ; do echo $(bjobs -w | grep PBSvsCyto.${ct} | wc -l) "jobs remaining"; sleep 5; done

# Check out if any jobs failed
cd ${LOG_ROOT}
rm -f ${LOG_ROOT}/temp.failedjobs.txt
for f in *out; do
	if ! grep -q "Successfully completed" ${f}; then
		echo "${f}" >> "${LOG_ROOT}/temp.failedjobs.txt"
	fi
done
echo "Total jobs: $(ls *.out | wc -l)  |    Failed jobs: $(wc -l temp.failedjobs.txt)"

# re-run failed jobs (53 failed. all due to memory limit)
CMD=${DIR}/step12_commands_rerun.txt
rm -f ${CMD}
while read f;do
	head -1 ${f} | sed 's/RUN: //g' >> ${CMD}
done < ${LOG_ROOT}/temp.failedjobs.txt
MAX_CONCURRENT=150
N=$(wc -l < ${CMD})
bsub -q long -W 2:00 -n 1 -R "span[hosts=1]" -R "rusage[mem=30000]" \
  -J "PBSvsCyto.${ct}[1-${N}]%${MAX_CONCURRENT}" \
  -o "${LOG_ROOT}/PBSvsCyto_%J_%I.out" \
  -e "${LOG_ROOT}/PBSvsCyto_%J_%I.err" \
  ${DIR}/run_task.sh ${CMD}
while [[ $(bjobs -w | grep PBSvsCyto.${ct} | wc -l) != 0 ]] ; do echo $(bjobs -w | grep PBSvsCyto.${ct} | wc -l) "jobs remaining"; sleep 5; done
cd ${LOG_ROOT}
rm -f ${LOG_ROOT}/temp.failedjobs.txt
for f in *49695*out; do
        if ! grep -q "Successfully completed" ${f}; then
                echo "${f}" >> "${LOG_ROOT}/temp.failedjobs.txt"
        fi
done
echo "Total jobs: $(ls *.out | wc -l)  |    Failed jobs: $(wc -l temp.failedjobs.txt)"



# After jobs run and QC'ed, summarize
for cytokine in IFNB IFNG TNF;do
	summary=${DIR}/coloc_susie/coloc_susie_summary_${cytokine}.txt
	rm -f ${summary}

	cd ${DIR}/coloc_susie/${cytokine}
	this_file=$(ls ${DIR}/coloc_susie/${cytokine} | grep -v "pdf" | grep IRF3 | head -1)
	awk '{OFS=FS="\t"}NR==1{print $0}' ${this_file} > ${summary}
	for f in coloc_susie*_summary.tsv; do
		awk '{OFS=FS="\t"}NR>1{print $0}' ${f} >> ${summary}
	done
	echo "${summary}"
done

Rscript ${DIR}/step12_coloc_susie_postRunCollection.R






