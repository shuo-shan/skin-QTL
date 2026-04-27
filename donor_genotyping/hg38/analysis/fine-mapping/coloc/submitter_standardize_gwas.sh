#!/bin/bash

dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc
cd $dir

# Create commands file
rm commands.txt
while read trait;do
	echo "bash ${dir}/standardize_gwas_for_coloc.sh ${trait}" >> commands.txt
done < traits.txt


while read trait;do
	cd ${dir}/${trait}
	f=$(ls *.h.tsv.gz)
	zcat ${f} | head -1 | tr '\t' '\n' > temp1
	zcat ${f} | head -2 | tail -1 | tr '\t' '\n' > temp2
	echo "------------------------------------------------------------"
	echo ${trait}
	paste temp1 temp2
	rm temp1 temp2
done < traits.txt



# Create a wrapper script
cat > ${dir}/run_task.sh << 'EOF'
#!/bin/bash
cmd=$(sed -n "${LSB_JOBINDEX}p" "$1")
echo "RUN: $cmd"
eval $cmd
EOF

chmod +x ${dir}/run_task.sh


# Submit jobs as an array
CMD=${dir}/commands.txt
N=$(wc -l < "$CMD")
MAX_CONCURRENT=400

mkdir -p ${dir}/log
log_dir=${dir}/log
bsub -q short -W 0:30 -n 1 -R "span[hosts=1]" -R "rusage[mem=1000]" \
  -J "parseByChr[1-${N}]%${MAX_CONCURRENT}" \
  -o "${dir}/log/parseByChr_%J_%I.out" \
  -e "${dir}/log/parseByChr_%J_%I.err" \
  ${dir}/run_task.sh $CMD


# Monitor jobs' progress
while [[ $(bjobs -w | grep parseByChr | wc -l) != 0 ]] ; do echo $(bjobs -w | grep parseByChr | wc -l) "jobs remaining"; sleep 5; done


