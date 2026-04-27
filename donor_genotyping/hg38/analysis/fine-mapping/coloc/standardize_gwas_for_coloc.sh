#!/bin/bash
# call Rscript to standardize GWAS stats data


dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc
trait=$1

echo "== ${trait} =="

cd ${dir}
trait=Melanomas_of_skin_dx_or_hx
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif Rscript standardize_gwas_for_coloc_${trait}.R
