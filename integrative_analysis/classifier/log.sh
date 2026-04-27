
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/classifier
cd ${dir}

##### 1. from CTCF-overlapping MEL reQTL (Pval1E-05), manually curate yes and no lists (based on shape of paired plots)

##### 2. copy the dictionary: 
ln -s /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/MELreQTL1E-05_overlapping_TFpeaks_with_modelingresults.txt

##### 3. copy the bed files of SNP and gene-level TSS
ln -s /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/MELreQTL1E-05_overlapping_ChIPpeaks.bed

