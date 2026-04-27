#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=300000]
#BSUB -q long
#BSUB -W 72:00
#BSUB -J plinkQC
#BSUB -e "./job_%J.err"
#BSUB -o "./job_%J.out"

### --------- in the cluster ----------- 
dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC
cd ${dir}

### first compile all necessary documents
### get test.vcf
module load bcftools
module load plink/1.90b6.27
# List all chromosomes in your VCF
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/merged_lowGPreset.vcf.gz

# create plink files of merged dataset
echo "creating plink files"; date
cd ${dir}/genodata
plink --vcf ${vcf} --make-bed --out data_nosex
echo "done creating plink files; annotating sex now"; date
# assign sex to plink .fam file
awk '{ $5 = 1; print }' data_nosex.fam > data_with_sex.fam
mv data_with_sex.fam data_nosex.fam
plink --bfile data_nosex --make-bed --out data
rm data_nosex*
plink --bfile data --pca --out data
plink --bfile data --sexcheck --out data
rm data_nosex*

# HapMap data and genotype PCA
cd ${dir}/genodata/qc

### ---------- Ancestry Estimation ---------------
qcdir=${dir}/genodata
refdir=${dir}/reference
name='data'
refname='HapMapIII_CGRCh38'
mkdir -p $qcdir/plink_log
cd ${qcdir}

# We will use an awk script to find A→T and C→G SNPs. As these SNPs are more difficult to align and only a subset of SNPs is required for the analysis, we will remove them from both the reference and study data set.
echo "step 1..."; date
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $qcdir/$name.bim  > \
    $qcdir/$name.ac_gt_snps

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $refdir/$refname.bim  > \
    $qcdir/$refname.ac_gt_snps
   
plink --bfile  $refdir/$refname \
      --exclude $qcdir/$refname.ac_gt_snps \
      --make-bed \
      --out $qcdir/$refname.no_ac_gt_snps
mv  $qcdir/$refname.no_ac_gt_snps.log $qcdir/plink_log/$refname.no_ac_gt_snps.log

plink --bfile  $qcdir/$name \
      --exclude $qcdir/$name.ac_gt_snps \
      --make-bed \
      --out $qcdir/$name.no_ac_gt_snps
mv  $qcdir/$name.no_ac_gt_snps.log $qcdir/plink_log/$name.no_ac_gt_snps.log

# Prune variants in linkage disequilibrium (LD) with an 𝑟2>0.2 in a 50kb window
echo "step2..."; date
cd $refdir; wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt \
     -O high-LD-regions-hg38-GRCh38.txt
highld=high-LD-regions-hg38-GRCh38.txt
plink --bfile  $qcdir/$name.no_ac_gt_snps \
      --exclude range  $refdir/$highld \
      --indep-pairwise 50 5 0.2 \
      --out $qcdir/$name.no_ac_gt_snps
mv  $qcdir/$name.no_ac_gt_snps.log $qcdir/plink_log

plink --bfile  $qcdir/$name.no_ac_gt_snps \
      --extract $qcdir/$name.no_ac_gt_snps.prune.in \
      --make-bed \
      --out $qcdir/$name.pruned
mv  $qcdir/$name.pruned.log $qcdir/plink_log/$name.pruned.log


# Use the list of pruned variants from the study sample to reduce the reference dataset to the size of the study samples:
echo "step3..."; date
plink --bfile  $refdir/$refname \
      --extract $qcdir/${name}.no_ac_gt_snps.prune.in \
      --make-bed \
      --out $qcdir/$refname.pruned
mv  $qcdir/$refname.pruned.log $qcdir/plink_log/$refname.pruned.log

# Check and correct chromosome mismatch
echo "step4..."; date
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim | \
    sed -n '/^[XY]/!p' > $qcdir/$refname.toUpdateChr

plink --bfile $qcdir/$refname.pruned \
      --update-chr $qcdir/$refname.toUpdateChr 1 2 \
      --make-bed \
      --out $qcdir/$refname.updateChr
mv $qcdir/$refname.updateChr.log $qcdir/plink_log/$refname.updateChr.log


# find variants with mis-matching chromosomal positions.
echo "step5..."; date
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
    $qcdir/${refname}.toUpdatePos

# Possible allele flips
echo "step6..."; date
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
    $qcdir/$refname.toFlip

# Upate positions and flip alleles
echo "step7..."; date
plink --bfile $qcdir/$refname.updateChr \
      --update-map $qcdir/$refname.toUpdatePos 1 2 \
      --flip $qcdir/$refname.toFlip \
      --make-bed \
      --out $qcdir/$refname.flipped
mv $qcdir/$refname.flipped.log $qcdir/plink_log/$refname.flipped.log

# remove mismatches
echo "step8..."; date
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.flipped.bim > \
    $qcdir/$refname.mismatch

plink --bfile $qcdir/$refname.flipped \
      --exclude $qcdir/$refname.mismatch \
      --make-bed \
      --out $qcdir/$refname.clean
mv $qcdir/$refname.clean.log $qcdir/plink_log/$refname.clean.log

# merge study genotype and reference data
echo "step9..."; date
# Step 1: Remove bad SNPs
plink --bfile $name.pruned \
      --exclude $name.merge.${refname}-merge.missnp \
      --make-bed \
      --out $name.pruned.nomiss

# Assign safe variant IDs
plink --bfile $name.pruned.nomiss \
      --set-missing-var-ids @:# \
      --make-bed \
      --out $name.pruned.nomiss.renamed

# Merge with HapMap
plink --bfile $name.pruned.nomiss.renamed \
      --bmerge $refname.clean \
      --make-bed \
      --out $qcdir/$name.merge.$refname

mv $qcdir/$name.merge.$refname.log $qcdir/plink_log

# PCA on the merged data
echo "step10..."; date
plink --bfile $qcdir/$name.merge.$refname \
      --pca \
      --out $qcdir/$name.$refname
mv $qcdir/$name.$refname.log $qcdir/plink_log


# this is for my QTL modeling: PCA on my pruned data
plink --bfile $qcdir/$name.pruned.nomiss.renamed \
      --pca 10 \
      --out $qcdir/$name.only
