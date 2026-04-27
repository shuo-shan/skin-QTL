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
traits_file=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/traits.txt
STEP8=${DIR}/step8_coloc_runRscript.sh
CMD_DIR=${DIR}/coloc/cmd_by_trait; mkdir -p ${CMD_DIR}
LOG_ROOT=${DIR}/logs/coloc;        mkdir -p ${LOG_ROOT}

gene_file=${DIR}/data/fdr05_genelist.txt
cat ${DIR}/eigenMT/results/*fdr05_genelist.txt | sort -u | sed '/^$/d' > ${gene_file}
#head -4 ${gene_file} > ${gene_file}.tmp && mv ${gene_file}.tmp ${gene_file}   # optional test

# create helper file
cat > ${DIR}/run_task.sh <<'EOF'
#!/bin/bash
set -euo pipefail
cmd=$(sed -n "${LSB_JOBINDEX}p" "$1"); echo "RUN: $cmd"; eval "$cmd"
EOF
chmod +x ${DIR}/run_task.sh

# Submit BSUB jobs
MAX_CONCURRENT=200

while read -r trait; do
  [[ -z "$trait" ]] && continue
  CMD=${CMD_DIR}/cmd_${trait}.txt
  awk -v ct="$ct" -v tr="$trait" -v step8="$STEP8" \
    '{print "bash "step8" "ct" "$0" "tr}' \
    ${gene_file} > ${CMD}

  N=$(wc -l < ${CMD})
  mkdir -p ${LOG_ROOT}/${trait}

  bsub -q short -W 0:30 -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" \
    -J "FRB_coloc_${trait}[1-${N}]%${MAX_CONCURRENT}" \
    -o "${LOG_ROOT}/${trait}/coloc_${trait}_%J_%I.out" \
    -e "${LOG_ROOT}/${trait}/coloc_${trait}_%J_%I.err" \
    ${DIR}/run_task.sh ${CMD}
done < ${traits_file}


# Monitor jobs' progress
while [[ $(bjobs -w | grep FRB_coloc | wc -l) != 0 ]] ; do echo $(bjobs -w | grep FRB_coloc | wc -l) "jobs remaining"; sleep 5; done


# ------- submit failed jobs for fibroblast only ------------ #
# because log file was not set-up properly in the first main run. this is a tailored solution for fibroblast
while read -r trait; do
	CMD=${CMD_DIR}/cmd_rerun_${trait}.txt

	# check which genes failed in plotting
	cd ${DIR}/coloc/${trait}/plots
	ls | sed 's/.locus.*//g' | sed 's/FRB_//g' | sort -u > genes_ok.txt
	comm -13 genes_ok.txt ${gene_file} > genes_failed.txt
	wc -l genes_failed.txt


	# make commands
	awk -v ct="$ct" -v tr="$trait" -v step8="$STEP8" \
	    '{print "bash "step8" "ct" "$0" "tr}' \
	    genes_failed.txt > ${CMD}
	rm genes_failed.txt genes_ok.txt
	head -1 ${CMD}


	# Submit BSUB jobs
      	N=$(wc -l < ${CMD})
	mkdir -p ${LOG_ROOT}/${trait}

	bsub -q short -W 2:00 -n 1 -R "span[hosts=1]" -R "rusage[mem=5000]" \
	  -J "FRB_coloc_${trait}[1-${N}]%50" \
	  -o "${LOG_ROOT}/${trait}/coloc_${trait}_%J_%I.out" \
	  -e "${LOG_ROOT}/${trait}/coloc_${trait}_%J_%I.err" \
	  ${DIR}/run_task.sh ${CMD}

done < ${traits_file}



# After jobs are done, check if all jobs are done
while read -r trait; do
        cd ${LOG_ROOT}/${trait}
        ls -lh | grep "Feb 16" | grep "out" | sed 's/.* //g' >> temp.rerunjobs.txt
        while read f; do
                grep "in cluster" ${f} | grep "Done" >> temp.rerun.donejobs.txt
        done < temp.rerunjobs.txt

        echo "Total jobs: $(wc -l temp.rerunjobs.txt)  |    Finished jobs: $(wc -l temp.rerun.donejobs.txt)    |   Trait: ${trait}"
        rm temp.rerun.donejobs.txt temp.rerunjobs.txt
done < ${traits_file}


# Check out why some jobs failed
while read -r trait; do
        cd ${LOG_ROOT}/${trait}
        mkdir -p ${LOG_ROOT}/${trait}
        rm -f ${LOG_ROOT}/${trait}/temp.failedjobs.txt
        for f in *out; do
                if ! grep -q "Done" "$f"; then
                        echo "${f}" >> "${LOG_ROOT}/${trait}/temp.failedjobs.txt"
                fi
        done
        echo "Total jobs: $(ls *.out | wc -l)  |    Failed jobs: $(wc -l temp.failedjobs.txt)    |   Trait: ${trait}"
done < ${traits_file}

# Most jobs failed due to memory issues (memory limit mostly, few run time limit), check reasons of other failed jobs
while read -r trait; do
        cd ${LOG_ROOT}/${trait}
        while read f;do
                if ! grep -q "memory usage limit" ${f}; then
                        echo ${f} >> temp.failedJobs.notCuzMemory.txt
                fi
        done < temp.failedjobs.txt
        cat temp.failedJobs.notCuzMemory.txt

done < ${traits_file}

# Create a list of failed command lines
while read -r trait; do
  cd ${LOG_ROOT}/${trait}
  if [ -f temp.failedjobs.txt ]; then
        CMD=${CMD_DIR}/rerun_commands_${trait}.txt
        rm -f ${CMD}

        # ------ Compile commands file ------ #
        while read -r f; do
              gene=$(grep "starting to perform coloc on" ${f} | sed 's/.* //g')
              echo "bash ${STEP8} ${ct} ${gene} ${trait}" >> ${CMD}
        done < temp.failedjobs.txt

        # ------ Submit BSUB jobs ------ #
        N=$(wc -l < ${CMD})
        echo ${N} ${CMD}
        mkdir -p ${LOG_ROOT}/rerun/${trait}

        bsub -q short -W 2:00 -n 1 -R "span[hosts=1]" -R "rusage[mem=5000]" \
          -J "coloc_${trait}[1-${N}]%50" \
          -o "${LOG_ROOT}/rerun/${trait}/coloc_${trait}_%J_%I.out" \
          -e "${LOG_ROOT}/rerun/${trait}/coloc_${trait}_%J_%I.err" \
          ${DIR}/run_task.sh ${CMD}
 fi
done < ${traits_file}

# ----------------------------------------------------------------------
# after jobs run, summarize
mkdir -p ${DIR}/coloc/summary
CMD=${CMD_DIR}/summarize_commands.txt
rm -f ${CMD}
QTLtype=eQTL
while read -r trait; do
        for cond in PBS IFNG IFNB TNF;do
                echo "bash ${DIR}/step8_coloc_summarize_results.sh ${ct} ${trait} ${cond} ${QTLtype}" >> ${CMD}
        done
done < ${traits_file}

QTLtype=reQTL
while read -r trait; do
        for cond in IFNG IFNB TNF;do
                echo "bash ${DIR}/step8_coloc_summarize_results.sh ${ct} ${trait} ${cond} ${QTLtype}" >> ${CMD}
        done
done < ${traits_file}

N=$(wc -l < ${CMD})
echo ${N} ${CMD}
mkdir -p ${DIR}/logs/coloc/summary

bsub -q short -W 0:30 -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" \
  -J "FRB.summary[1-${N}]%100" \
  -o "${DIR}/logs/coloc/summary/coloc_summary_%J_%I.out" \
  -e "${DIR}/logs/coloc/summary/coloc_summary_%J_%I.err" \
  ${DIR}/run_task.sh ${CMD}



# ------------------------------------------------------------------------
# Join coloc summary by QTL and GWAS statistics
script=${DIR}/step8_coloc_postRunSummary_runRscript.sh
while read -r trait; do
  CMD=${CMD_DIR}/cmd_join_coloc_summary_by_stats_commands_${trait}.txt
  rm -f ${CMD}
  awk -v ct="$ct" -v tr="$trait" -v script="$script" \
    '{print "bash "script" "ct" "$0" "tr}' \
    ${gene_file} > ${CMD}

  N=$(wc -l < ${CMD})
  mkdir -p ${LOG_ROOT}/join/ ${LOG_ROOT}/join/${trait}
  echo ${trait} ${N}

  bsub -q short -W 2:00 -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" \
    -J "FRB_coloc_${trait}[1-${N}]%200" \
    -o "${LOG_ROOT}/coloc/join/${trait}/coloc_${trait}_%J_%I.out" \
    -e "${LOG_ROOT}/coloc/join/${trait}/coloc_${trait}_%J_%I.err" \
    ${DIR}/run_task.sh ${CMD}
done < ${traits_file}


# Check out if any jobs failed
while read -r trait; do
        cd ${LOG_ROOT}/coloc/join/${trait}
        rm -f ${LOG_ROOT}/coloc/join/${trait}/temp.failedjobs.txt
	rm -f ${LOG_ROOT}/coloc/join/${trait}/temp.maxmemory.txt
        for f in *out; do
		grep "Max Memory" ${f} >> "${LOG_ROOT}/coloc/join/${trait}/temp.maxmemory.txt"
		if ! grep -q "Done" "$f"; then
         	       echo "${f}" >> "${LOG_ROOT}/coloc/join/${trait}/temp.failedjobs.txt"
	        fi
        done
        echo "Total jobs: $(ls *.out | wc -l)  |    Failed jobs: $(wc -l temp.failedjobs.txt)    |   Trait: ${trait}"
done < ${traits_file}

# Most jobs failed due to memory issues (memory limit mostly, few run time limit), check reasons of other failed jobs
while read -r trait; do
	cd ${LOG_ROOT}/coloc/join/${trait}
        while read f;do
                if ! grep -q "memory usage limit" ${f}; then
                        echo ${f} >> temp.failedJobs.notCuzMemory.txt
                fi
        done < temp.failedjobs.txt
        cat temp.failedJobs.notCuzMemory.txt

done < ${traits_file}

# Create a list of failed command lines and re-run
while read -r trait; do
  cd ${LOG_ROOT}/coloc/join/${trait}
  if [ -f temp.failedjobs.txt ]; then
	  mkdir -p ${LOG_ROOT}/rerun/join/${trait}
        CMD=${CMD_DIR}/rerun_join_commands_${trait}.txt
        rm -f ${CMD}

        # ------ Compile commands file ------ #
        while read -r f; do
	      script=${DIR}/step8_coloc_postRunSummary_runRscript.sh
              gene=$(grep "starting to perform coloc on" ${f} | sed 's/.* //g')
              echo "bash ${script} ${ct} ${gene} ${trait}" >> ${CMD}
        done < temp.failedjobs.txt

        # ------ Submit BSUB jobs ------ #
        N=$(wc -l < ${CMD})
        echo ${N} ${CMD}
        mkdir -p ${LOG_ROOT}/rerun/${trait}

        bsub -q short -W 2:00 -n 1 -R "span[hosts=1]" -R "rusage[mem=5000]" \
          -J "coloc_${trait}[1-${N}]%100" \
          -o "${LOG_ROOT}/rerun/join/${trait}/coloc_${trait}_%J_%I.out" \
          -e "${LOG_ROOT}/rerun/join/${trait}/coloc_${trait}_%J_%I.err" \
          ${DIR}/run_task.sh ${CMD}
 fi
done < ${traits_file}

#  ------------------------------------------------------------------------
# summarize joined coloc summary and statistics of QTL and GWAS
mkdir -p ${DIR}/coloc/summary
CMD=${CMD_DIR}/summarize_commands.txt
rm -f ${CMD}
QTLtype=eQTL
while read -r trait; do
        for cond in PBS IFNG IFNB TNF;do
                echo "bash ${DIR}/step8_coloc_summarize_joined_results.sh ${ct} ${trait} ${cond} ${QTLtype}" >> ${CMD}
        done
done < ${traits_file}

QTLtype=reQTL
while read -r trait; do
        for cond in IFNG IFNB TNF;do
                echo "bash ${DIR}/step8_coloc_summarize_joined_results.sh ${ct} ${trait} ${cond} ${QTLtype}" >> ${CMD}
        done
done < ${traits_file}

N=$(wc -l < ${CMD})
echo ${N} ${CMD}
mkdir -p ${DIR}/logs/coloc/summary

bsub -q short -W 0:30 -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" \
  -J "${ct}_summary[1-${N}]%100" \
  -o "${DIR}/logs/coloc/summary/coloc_summary_%J_%I.out" \
  -e "${DIR}/logs/coloc/summary/coloc_summary_%J_%I.err" \
  ${DIR}/run_task.sh ${CMD}
