#!/bin/bash

dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/manhattan_plot
#mkdir ${dir}/log
# get SNPs with a calculated p value for PBS-eQTL featureselected model
cd ${dir}
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all.txt
cat ${f} | awk '{if ($15 != ".") print $0}' > temp_result.txt
cat temp_result.txt | cut -f1 | sort | uniq | shuf > temp_snps.txt

# for each SNP, find the gene with the strongest association and store.
# 706390 temp_snps.txt
# split files to 1500 snps per file, total of 471 files.
jobname=cs
mkdir split_snps; cd split_snps; split -l 35000 ${dir}/temp_snps.txt splitSnps_
for this_file in splitSnps_*;do
	rm ${dir}/commands.txt ${dir}/shuffled_commands.txt ${dir}/commands.joined.txt
	cat ${dir}/split_snps/${this_file} | awk '{OFS=" "}{print "bash /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/manhattan_plot/script_step2.sh",$0,"/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/manhattan_plot"}' > ${dir}/commands.txt
	sleep 3
	shuf ${dir}/commands.txt > ${dir}/shuffled_commands.txt
	bash /pi/manuel.garber-umw/sshan/scripts/function_collapse_commands.sh ${dir} ${dir}/shuffled_commands.txt ${dir}/commands.joined.txt 100
	sleep 3
	while read c; do
		echo ${c} | bsub -W 05:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=400] -q long -e "${dir}/log/step2_%J%I.err" -o "${dir}/log/step2_%J%I.out"
	done < ${dir}/commands.joined.txt

	while [[ $(bjobs | grep ${jobname} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${jobname} | wc -l) "jobs remaining"; date;sleep 5; done

	# clean up
	mkdir ${dir}/log/${this_file}; mv ${dir}/log/step2_* ${dir}/log/${this_file}
done

