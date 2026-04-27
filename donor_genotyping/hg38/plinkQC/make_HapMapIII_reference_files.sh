#!/bin/bash
# script adapted from: https://meyer-lab-cshl.github.io/plinkQC/articles/HapMap.html

### get test.vcf
module load plink/1.90b6.27

### ---------  make HapMap reference plink files  ---------  
dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC
refdir=${dir}/reference
mkdir -p ${refdir}
mkdir -p ${dir}/genodata/qc/plink_log

cd ${refdir}

ftp=ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/
prefix=hapmap3_r2_b36_fwd.consensus.qc.poly

wget $ftp/$prefix.map.bz2
bunzip2 $prefix.map.bz2

wget $ftp/$prefix.ped.bz2
bunzip2 $prefix.per.bz2

wget $ftp/relationships_w_pops_121708.txt

plink --file $refdir/$prefix \
      --make-bed \
      --out $refdir/HapMapIII_NCBI36
mv $refdir/HapMapIII_NCBI36.log $refdir/log

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz
gunzip hg18ToHg38.over.chain.gz
awk '{print "chr" $1, $4 -1, $4, $2 }' $refdir/HapMapIII_NCBI36.bim | \
    sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > \
    $refdir/HapMapIII_NCBI36.tolift
conda create -n liftover -c bioconda ucsc-liftover
conda activate liftover
liftOver $refdir/HapMapIII_NCBI36.tolift $refdir/hg18ToHg38.over.chain \
    $refdir/HapMapIII_CGRCh38 $refdir/HapMapIII_NCBI36.unMapped
# extract i) the variants that were mappable from the old to the new genome and ii) their updated positions
# ectract mapped variants
awk '{print $4}' $refdir/HapMapIII_CGRCh38 > $refdir/HapMapIII_CGRCh38.snps
# ectract updated positions
awk '{print $4, $3}' $refdir/HapMapIII_CGRCh38 > $refdir/HapMapIII_CGRCh38.pos
# use PLINK to extract the mappable variants from the old build and update their position
plink --bfile $refdir/HapMapIII_NCBI36 \
    --extract $refdir/HapMapIII_CGRCh38.snps \
    --update-map $refdir/HapMapIII_CGRCh38.pos \
    --make-bed \
    --out $refdir/HapMapIII_CGRCh38
mv $refdir/HapMapIII_CGRCh38.log $refdir/log

