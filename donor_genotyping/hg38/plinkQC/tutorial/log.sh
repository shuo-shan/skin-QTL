
### --------- in the cluster ----------- 
dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC
cd ${dir}
mkdir -p genodata
mkdir -p genodata/qc

### first compile all necessary documents
### get test.vcf
module load bcftools
module load plink/1.90b6.27
# List all chromosomes in your VCF
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/merged_lowGPreset.vcf.gz
bcftools view -h ${vcf} | grep "^##contig"
# Loop through each chromosome and sample 10,000 variants
outdir=${dir}/sampled_by_chr
mkdir -p ${outdir}

for chr in chr{1..22} chrX chrY; do
  echo "Processing $chr"

  bcftools view -r $chr $vcf \
    | bcftools view -H \
    | shuf -n 10000 \
    | cut -f1,2 \
    > ${outdir}/${chr}_random10000.txt
done

cat ${outdir}/chr*_random10000.txt > ${outdir}/sampled_positions.txt
bcftools view -T ${outdir}/sampled_positions.txt ${vcf} -Oz -o sampled_merged.vcf.gz
bcftools index sampled_merged.vcf.gz

# Create PCA of merged dataset:
cd ${dir}/genodata
plink --vcf ${dir}/sampled_merged.vcf.gz --make-bed --out test_nosex # generate .bed .bim .fam

# assign sex to plink .fam file
cd ${dir}/genodata
awk '{ $5 = 1; print }' test_nosex.fam > test_with_sex.fam
mv test_with_sex.fam test_nosex.fam
plink --bfile test_nosex --make-bed --out test
rm test_nosex*
plink --bfile test --pca --out test
plink --bfile test --freq --out test
plink --bfile test --sexcheck --out test

rm test_nosex* test_with_sex*

# HapMap data and genotype PCA
cd ${dir}/genodata/qc
wget https://github.com/meyer-lab-cshl/plinkQC/blob/master/tests/testthat/HapMap_ID2Pop.txt
wget https://github.com/meyer-lab-cshl/plinkQC/blob/master/tests/testthat/HapMap_PopColors.txt
wget https://github.com/meyer-lab-cshl/plinkQC/blob/master/tests/testthat/data.HapMapIII.eigenvec


### --------- Rscript  ----------- 
module load r/4.2.2
R
# Rscript
library(plinkQC)
indir="/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata"
qcdir="/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata/qc"
name <- "test" # # Because your files are test.bed, test.bim, test.fam
path2plink <- "/share/pkg/plink/1.90b6.27/plink"

fail_individuals <- perIndividualQC(
  indir = indir,
  qcdir = qcdir,
  name = name,
  path2plink = path2plink,
  refSamplesFile = file.path(qcdir, "HapMap_ID2Pop.txt"),
  refColorsFile = file.path(qcdir, "HapMap_PopColors.txt"),
  prefixMergedDataset = file.path(qcdir, "data.HapMapIII"),  # only if you downloaded this .bed/.bim/.fam set (not just eigenvec)
  do.run_check_ancestry = TRUE,
  interactive = TRUE,
  verbose = TRUE
)

# overview individual level QC
overview_individuals <- overviewPerIndividualQC(fail_individuals,
                                                interactive=TRUE)
### ---------  make HapMap reference plink files  ---------  
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



### ---------- Ancestry Estimation ---------------

qcdir==${dir}/genodata
refdir=${dir}/reference
name='data'
refname='HapMapIII'

mkdir -r $qcdir/plink_log






