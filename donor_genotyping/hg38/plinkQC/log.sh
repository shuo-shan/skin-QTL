
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

# Create PCA of merged dataset:
cd ${dir}/genodata
plink --vcf ${vcf} --make-bed --out data_nosex
# assign sex to plink .fam file
awk '{ $5 = 1; print }' data_nosex.fam > data_with_sex.fam
mv data_with_sex.fam data_nosex.fam
plink --bfile data_nosex --make-bed --out data
rm test_nosex*
plink --bfile data --pca --out data
plink --bfile data --sexcheck --out data
rm data_nosex*

# HapMap data and genotype PCA
cd ${dir}/genodata/qc
wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/HapMap_ID2Pop.txt -O HapMap_ID2Pop.txt
wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/HapMap_PopColors.txt -O HapMap_PopColors.txt
wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/data.HapMapIII.eigenvec -O data.HapMapIII.eigenvec


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




