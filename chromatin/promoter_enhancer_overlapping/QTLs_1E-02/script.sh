#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=18000]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e “./%J%I.err" 
#BSUB -o “./%J%I.out"

module load bedtools
inFile=$1
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/promoter_enhancer_overlapping/QTLs_1E-02

while read line;do
	tag=$(echo ${line} | tr ' ' '_')
        echo ${line} | tr ' ' '\t' > this_QTL_${tag}.bed
        this_n=$(bedtools intersect -loj -a this_QTL_${tag}.bed -b CRE_dynamics_3cts.bed | grep 'promoter\|enhancer' | wc -l)

        if [ ${this_n} == 0 ];then
                bedtools intersect -loj -a this_QTL_${tag}.bed -b CRE_dynamics_3cts.bed | head -1 > result_${tag}.txt 
        else
                bedtools intersect -loj -a this_QTL_${tag}.bed -b CRE_dynamics_3cts.bed | grep 'promoter\|enhancer' | head -1 > result_${tag}.txt
        fi      

	output_file="${dir}/QTL_overlapping_CRE.txt"
	lock_file="${dir}/QTL_overlapping_CRE.lock"
	exec 9>"$lock_file" # opens the lock file for writing and associates it with file descriptor 9. lock file is created if it doesn't already exist.
	flock 9 # apply an advisory lock on this file
	cat result_${tag}.txt >> ${output_file}
	# release the lock
	flock -u 9

	# cleanup
	rm result_${tag}.txt this_QTL_${tag}.bed
done < ${inFile}
