# log of commands used to trim reads. reads were previously trimmed 3' end from 150bp to 100bp to prevent mapping to polyA. reads were also barcode modified to cater to ESAT requirements.

Dir=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/fq_read_length_comparison
cd ${Dir}
mkdir -p fq_40bp log

# generate commands

# bsub command
while read c; do
   echo ${c} | bsub -W 1:00 -n 1 -R "span[hosts=1]" -R rusage[mem=1020] -q long -e "./log/trim_%J.err" -o "./log/trim_%J.out"
done < commands.txt
