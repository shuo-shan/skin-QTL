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

# create batch job commands
echo "creating commands for QTL genes"
cd ${DIR}/eigenMT/results/ 
cat *fdr05_genelist.txt | sort -u > fdr05_genelist.txt

awk '{print "bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/step7_susie_runRscript.sh",$0}' ${DIR}/eigenMT/results/fdr05_genelist.txt > ${DIR}/susie/commands_susie.txt

# submit BSUB jobs
CMD_FILE=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/susie/commands_susie.txt
#CMD_FILE=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/commands_failed_files.txt
LOG_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs

N=$(wc -l < "$CMD_FILE")

bsub -J "susie[1-${N}]%20" \
  -q short \
  -n 1 \
  -R "rusage[mem=20000]" \
  -M 20000 \
  -W 0:50 \
  -oo "${LOG_DIR}/susie_%J_%I.out" \
  -eo "${LOG_DIR}/susie_%J_%I.err" \
  "set -euo pipefail; cmd=\$(sed -n \"\${LSB_JOBINDEX}p\" $CMD_FILE); echo \"RUN: \$cmd\"; eval \"\$cmd\""

cd ${DIR}
