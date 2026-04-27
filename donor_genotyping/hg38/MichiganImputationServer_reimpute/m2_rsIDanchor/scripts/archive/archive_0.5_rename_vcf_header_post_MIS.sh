#!/bin/bash
#BSUB -J reheader_vcf
#BSUB -R "rusage[mem=64000]"
#BSUB -o reheader_vcfs_%I.out
#BSUB -e reheader_vcfs_%I.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1

# Load necessary modules
module load plink2/alpha6.1amd
module load bcftools
module load htslib
module load plink

DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/MIS_results
cd ${DIR}
INPUT_DIRS=("n101_1000G" "n101_HLA" "n101_HRC")


### ======== rename header in all vcf.gz files to be consistent throughout batches ====== ###
### n101_1000G
cd ${INPUT_DIRS[0]} 
bcftools query -l chr6.dose.vcf.gz | \
awk '{
    old = $0
    new = old
    sub(/^0_skineQTL-/, "", new)
    sub(/^0_/, "", new)
    sub(/^F0/, "F", new)
    print old "\t" new
}' > rename_samples.txt

for chr in {1..22}; do
	echo ${chr}
	bcftools reheader -s rename_samples.txt \
		-o chr${chr}.dose.tmp.vcf.gz \
		chr${chr}.dose.vcf.gz
	mv chr${chr}.dose.tmp.vcf.gz chr${chr}.dose.vcf.gz
	tabix -f -p vcf chr${chr}.dose.vcf.gz

	bcftools reheader -s rename_samples.txt \
                -o chr${chr}.empiricalDose.tmp.vcf.gz \
		chr${chr}.empiricalDose.vcf.gz
	mv chr${chr}.empiricalDose.tmp.vcf.gz chr${chr}.empiricalDose.vcf.gz
	tabix -f -p vcf chr${chr}.empiricalDose.vcf.gz
done

### n101_HRC
cd ${DIR}/${INPUT_DIRS[2]}
bcftools query -l chr6.dose.vcf.gz | \
awk '{
    old = $0
    new = old
    sub(/^0_skineQTL-/, "", new)
    sub(/^0_/, "", new)
    sub(/^F0/, "F", new)
    print old "\t" new
}' > rename_samples.txt

for chr in {1..22}; do
        echo ${chr}
        bcftools reheader -s rename_samples.txt \
                -o chr${chr}.dose.tmp.vcf.gz \
                chr${chr}.dose.vcf.gz
        mv chr${chr}.dose.tmp.vcf.gz chr${chr}.dose.vcf.gz
        tabix -f -p vcf chr${chr}.dose.vcf.gz

        bcftools reheader -s rename_samples.txt \
                -o chr${chr}.empiricalDose.tmp.vcf.gz \
                chr${chr}.empiricalDose.vcf.gz
        mv chr${chr}.empiricalDose.tmp.vcf.gz chr${chr}.empiricalDose.vcf.gz
        tabix -f -p vcf chr${chr}.empiricalDose.vcf.gz
done


### n101_HLA
cd ${DIR}/${INPUT_DIRS[1]}
bcftools query -l chr6.dose.vcf.gz | \
awk '{
    old = $0
    new = old
    sub(/^skineQTL-/, "", new)
    sub(/^F0/, "F", new)
    print old "\t" new
}' > rename_samples.txt

chr=6
bcftools reheader -s rename_samples.txt \
        -o chr${chr}.dose.tmp.vcf.gz \
        chr${chr}.dose.vcf.gz   
mv chr${chr}.dose.tmp.vcf.gz chr${chr}.dose.vcf.gz
tabix -f -p vcf chr${chr}.dose.vcf.gz











