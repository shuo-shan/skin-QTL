while read c; do
  echo ${c} | bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=20400] -R "select[rh=6]" -q long -e "./%J%I.err" -o "./%J%I.out"
done < commands.txt
