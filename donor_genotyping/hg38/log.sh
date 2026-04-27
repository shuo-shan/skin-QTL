
#
dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38
cd ${dir}

# look for all vcf files recursively
# I removed the VB and CB samples. only kept skin-eQTL project related samples in files.txt
find "$(pwd)" -type f -name "*.vcf.gz" -exec bash -c '
for f; do
  full=$(realpath "$f")
  dir=$(dirname "$f")
  base=$(basename "$f" .vcf.gz)
  echo -e "$full\t$dir\t$base"
done
' _ {} + > files.txt

# for each file, set any genotype of MAX(GP)<0.90 to be unknown, ./.
bash reset_lowMaxGP_to_unknown.sh ${c}
# arg1: full input vcf path
# arg2: fname prefix
# arg3: dir
cat files.txt | awk '{print "bash reset_lowMaxGP_to_unknown.sh",$0}' > commands_resetGT.txt

while read c; do
  echo ${c} | bsub -J resetGT -W 3:00 -n 1 -R "span[hosts=1]" -R rusage[mem=3000] -q long -e "./log/resetGT.%J.err" -o "./log/resetGT.%J.out"
done < commands_resetGT.txt


# find all lowGPresetGT.vcf.gz files
find "$(pwd)" -type f -name "*lowGPresetGT.vcf.gz" > vcf_list.txt
c="module load bcftools; dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38; bcftools merge -l vcf_list.txt -Oz -o ${dir}/merged_lowGPreset.vcf.gz; bcftools index ${dir}/merged_lowGPreset.vcf.gz"
echo ${c} | bsub -J merge -W 24:00 -n 1 -R "span[hosts=1]" -R rusage[mem=300GB] -q long -e "./log/merge.%J.err" -o "./log/merge.%J.out"


### ---------- make test file
bcftools view merged_lowGPreset.vcf.gz | head -10000 > test.vcf
bcftools view test.vcf -Oz -o test.vcf.gz
bcftools index test.vcf.gz

### ---------- plinkQC -------------
module load plink/1.90b6.27
conda create -n plinkQC -c bioconda -c conda-forge plink r-base
### per individual QC


### per variant QC



### generate cleaned data
