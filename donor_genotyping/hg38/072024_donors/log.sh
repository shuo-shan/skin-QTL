# record of commands used to download, process, and analysis of 07/2024 donor genotyping data
# shuo.shan@umassmed.edu, 11/2024

dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors

# download imputed vcf and metadata (QC, coverage, ancestry)
cd ${dir}
gencove download ./download/ --project-id f52abc44-3def-4583-b878-0379e7f0c439 --file-types impute-vcf
gencove download ./download/ --project-id f52abc44-3def-4583-b878-0379e7f0c439 --file-types ancestry-json,cnv-png,impute-csi,qc

# organize into respective folders
cd ${dir}
find . -type f -name "*vcf.gz*" > file_list.txt
awk '{print "mv "$0" ../organized"}' file_list.txt > command.txt; bash command.txt

find . -type f -name "*ancestry-json.json" > file_list.txt
awk '{print "mv "$0" ../ancestry"}' file_list.txt > command.txt; bash command.txt

find . -type f -name "*_qc.json" > file_list.txt
awk '{print "mv "$0" ../QC"}' file_list.txt > command.txt; bash command.txt
python ${dir}/QC_json_to_table.py
# coverage is plotted by ${dir}/overview_genotyping_coverage.R and saved as ${dir}/coverage_summary.png

find . -type f -name "*_cnv-png.png" > file_list.txt
awk '{print "mv "$0" ../CNV"}' file_list.txt > command.txt; bash command.txt


# modify so that all lowquality SNPs genotype is changed to ./. instead of 0/0
##FILTER=<ID=LOWCONF,Description="Set if not true: MAX(FORMAT/GP)>0.9">
cd /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors/organized
ls | grep -v gz.csi > file_list.txt
rm commands.txt
while read f;do
	echo "cd /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors/organized; module load bcftools; bcftools +setGT ${f} -- -t q -i 'FILTER=\"LOWCONF\"' -n './.' | bcftools view -Oz -o modified_${f}; bcftools index modified_${f}; echo ${f}" >> commands.txt
done < file_list.txt

while read c; do
  echo ${c} | bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10200] -q long -e "./%J%I.err" -o "./%J%I.out" -J modify
done < commands.txt

mv modified_* ${dir}/organized_lowSNPqualityGTModified
# once done merging, these files are deleted to save space, regenerating takes about 10 minutes using the above code


# fetch file ID and donor ID
cd ${dir}/organized
for f in *.vcf.gz;do
	donor=$(bcftools view -h ${f} | tail -1 | cut -f10)
	echo -e "${f}\t${donor}" >> ${dir}/file_donor_lookup.txt
done


# merge all samples
cd ${dir}/organized_lowSNPqualityGTModified
echo "cd ${dir}/organized_lowSNPqualityGTModified; ls *.vcf.gz > file_list.txt; bcftools merge -l file_list.txt -Oz -o ${dir}/merged/merged.vcf.gz; bcftools index merged.vcf.gz" > command.txt
while read c; do
  echo ${c} | bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=100200] -q long -e "./%J%I.err" -o "./%J%I.out" -J merge
done < command.txt


















