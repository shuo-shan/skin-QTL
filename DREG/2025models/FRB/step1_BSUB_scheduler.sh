#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=200]
#BSUB -q long
#BSUB -W 12:00
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/step1_scheduler_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/step1_scheduler_%J_%I.err"


echo "sleeping now"
sleep 1000
echo "slept for 1000 seconds, 4000 to go"
sleep 1000
echo "slept for 2000 seconds, 3000 to go"
sleep 1000
echo "slept for 3000 seconds, 2000 to go"
sleep 1000
echo "slept for 4000 seconds, 1000 to go"
sleep 1000
echo "slept for 5000 seconds, submitting jobs now"
cd /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB
bsub < step1_submit_chunks.sh
