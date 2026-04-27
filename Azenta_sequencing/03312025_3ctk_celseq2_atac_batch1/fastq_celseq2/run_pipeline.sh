while read c; do
   echo ${c} | bsub -W 12:00 -n 1 -R "span[hosts=1]" -R rusage[mem=1020] -q long -e "./log/pipeline_%J.err" -o "./log/pipeline_%J.out"
done < commands.txt


while read c; do
   echo ${c} | bsub -W 12:00 -n 1 -R "span[hosts=1]" -R rusage[mem=1020] -q long -e "./log/pipeline_demux_%J.err" -o "./log/pipeline_demux_%J.out"
done < commands.txt


while read c; do
   echo ${c} | bsub -W 12:00 -n 1 -R "span[hosts=1]" -R rusage[mem=1020] -q long -e "./log/pipeline_bcmod_%J.err" -o "./log/pipeline_bcmod_%J.out"
done < commands.txt
