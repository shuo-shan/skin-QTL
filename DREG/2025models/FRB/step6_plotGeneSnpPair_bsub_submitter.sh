#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=1000]
#BSUB -q long
#BSUB -W 08:00
#BSUB -J job_submitter
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/job_submitter_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/job_submitter_%J.err"
# batch submit jobs to plot gene:snp pairs

# reduces heavy I/O thundering-herd effect in parallel jobs
sleep $((RANDOM % 20))

# set-up
ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}

# specify fdr05 gene table
f=FRB_TNF_reQTL_gene_fdr05_table.txt
basename=$(echo ${f} | sed 's/\_table\.txt//g')
mkdir -p ${DIR}/plots ${DIR}/plots/temp_output
mkdir -p ${DIR}/data

## ---------- create gene chunk lookup table ---------- #
#cd ${DIR}/chunks
#echo -e "chunk\tgene" > ${DIR}/data/gene_chunk_dict.txt
#for f in gene_chunk_*;do
#       id=$(echo ${f} | cut -d"_" -f3)
#       echo ${f} " ID is: " ${id}
#       awk -v id=${id} '{OFS="\t"}{print id,$0}' ${f} >> ${DIR}/data/gene_chunk_dict.txt
#done
#echo "wrote gene chunk dict to ${DIR}/data/gene_chunk_dict.txt"

# create batch job commands
echo "creating commands for ${f}"
awk 'NR>1' ${DIR}/eigenMT/results/${f} | cut -f4,9 | tr "\t" " " > ${DIR}/plots/pairs_for_${f}
awk '{print "bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/step6_plotGeneSnpPair.sh",$0}' ${DIR}/plots/pairs_for_${f} > ${DIR}/plots/commands_for_${f}
#head -20 ${DIR}/plots/commands_for_${f} > ${DIR}/plots/top20.commands_for_${f}

# submit jobs as an array, concurrency = 300 jobs
CMD=${DIR}/plots/commands_for_${f}
N=$(wc -l < "$CMD")

bsub -q short -W 0:30 -n 1 -R "span[hosts=1]" -R "rusage[mem=500]" \
  -J "plot_${ct}[1-${N}]%300" \
  -o "${DIR}/logs/plot_%J_%I_${f}.out" \
  -e "${DIR}/logs/plot_%J_%I_${f}.err" \
  bash -lc '
    cmd=$(sed -n "${LSB_JOBINDEX}p" "'"$CMD"'")
    echo "RUN: $cmd"
    eval "$cmd"
  '

# monitor jobs' progress
while [[ $(bjobs -w | grep plot_${ct} | wc -l) != 0 ]] ; do echo $(bjobs -w | grep plot_${ct} | wc -l) "jobs remaining"; sleep 5; done

# compile all PDF into one big PDF
cd ${DIR}/plots/temp_output
mkdir -p ${basename}
mv *.pdf ${basename}

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step6_combine_all_plots.R ${DIR}/plots/temp_output/${basename}

cd ${DIR}/plots
rm ${DIR}/plots/pairs_for_${f} ${CMD}
#rm -r ${DIR}/plots/{A..Z}
rm -r ${DIR}/plots/temp_output
