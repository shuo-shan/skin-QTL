#!/bin/bash
# download GWAS summary statistics for each trait

dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/

# Record number of cases and controls per trait
echo -e "trait\tcases\tcontrols" > ${dir}/summary_all_traits.txt
echo -e "vitiligo_verma2024\t1263\t449492" >> ${dir}/summary_all_traits.txt
echo -e "vitiligo_jin2016\t2853\t37405" >> ${dir}/summary_all_traits.txt
echo -e "psoriasis\t36466\t458078" >> ${dir}/summary_all_traits.txt
echo -e "cutaneous_lupus_erythematosus\t752\t450411" >> ${dir}/summary_all_traits.txt



# Vitiligo
# 2024 Verma study (1,263 European ancestry cases, 449,492 European ancestry controls)
trait=vitiligo_verma2024
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476174/harmonised/GCST90476174.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476174/harmonised/GCST90476174.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476174/harmonised/GCST90476174.h.tsv.gz-meta.yaml

# Vitiligo: 2016 Jin study (2,853 European ancestry cases and 37,405 European ancestry controls)
trait=vitiligo_jin2016
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
# already downloaded before, GCST004785

# Psoriasis
# 2025 largest meta analysis (36,466 European ancestry cases, 458,078 European ancestry controls)
trait=psoriasis
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
# harnomized GWAS
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90472001-GCST90473000/GCST90472771/harmonised/GCST90472771.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90472001-GCST90473000/GCST90472771/harmonised/GCST90472771.h.tsv.gz-meta.yaml

# Cutaneous lupus erythematosus
# 2024 Verma study (752 European ancestry cases, 450,411 European ancestry controls)
trait=cutaneous_lupus_erythematosus
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476182/harmonised/GCST90476182.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476182/harmonised/GCST90476182.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476182/harmonised/GCST90476182.h.tsv.gz-meta.yaml

# Systemic lupus erythematosus
# 2024 Verma study (1,013 European ancestry cases, 449,940 European ancestry controls)
trait=systemic_lupus_erythematosus
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476183/harmonised/GCST90476183.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476183/harmonised/GCST90476183.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476183/harmonised/GCST90476183.h.tsv.gz-meta.yaml


# basal cell carcinoma
# 2024 Verma veteran study (31142 cases, 404406 controls, European)
trait=basal_cell_carcinoma
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
# harmonized GWAS stats, index, metadata
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475582/harmonised/GCST90475582.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475582/harmonised/GCST90475582.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475582/harmonised/GCST90475582.h.tsv.gz-meta.yaml


# squamous cell carcinoma
# 2024 veteran study (19217 cases, 419674 controls, European)
trait=squamous_cell_carcinoma
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
# harmonized GWAS stats, index, metadata
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475583/harmonised/GCST90475583.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475583/harmonised/GCST90475583.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475583/harmonised/GCST90475583.h.tsv.gz-meta.yaml


# rheumatoid arthritis
# 2024 Verma study (25,533 European ancestry cases, 290,135 European ancestry controls)
trait=rheumatoid_arthritis
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479433/harmonised/GCST90479433.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479433/harmonised/GCST90479433.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479433/harmonised/GCST90479433.h.tsv.gz-meta.yaml


# alopecia_areata
# 2024 Verma Science Veteran study (600 European ancestry cases, 450,317 European ancestry controls)
trait=alopecia_areata
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90478001-GCST90479000/GCST90478824/harmonised/GCST90478824.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90478001-GCST90479000/GCST90478824/harmonised/GCST90478824.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90478001-GCST90479000/GCST90478824/harmonised/GCST90478824.h.tsv.gz-meta.yaml


# melanoma
# 2024 Verma Science Veteran study 
trait=Melanomas_of_skin_dx_or_hx
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475577/harmonised/GCST90475577.h.tsv.gz


# atopic_dermatitis
# 2025 Nat Comm Olivia M Atopic dermatitis study (42,963 European ancestry cases, 408,472 European ancestry controls)
trait=atopic_dermatitis
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90503001-GCST90504000/GCST90503109/harmonised/GCST90503109.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90503001-GCST90504000/GCST90503109/harmonised/GCST90503109.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90503001-GCST90504000/GCST90503109/harmonised/GCST90503109.h.tsv.gz-meta.yaml


# Crohn's disease
# 2024 Verma Science Veteran study (2,256 European ancestry cases, 313,412 European ancestry controls)
trait=crohns_disease
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475318/harmonised/GCST90475318.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475318/harmonised/GCST90475318.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475318/harmonised/GCST90475318.h.tsv.gz-meta.yaml


# Sunburn
# 2018 Loh Nat Genet study (350,232 European ancestry individuals), hg38
trait=sunburn
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90029001-GCST90030000/GCST90029034/harmonised/29892013-GCST90029034-EFO_0003958.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90029001-GCST90030000/GCST90029034/harmonised/29892013-GCST90029034-EFO_0003958.h.tsv.gz-meta.yaml


# Skin pigmentation 
# 2025 UKbiobank (415,018 European ancestry individuals)
trait=skin_pigmentation
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90691001-GCST90692000/GCST90691754/harmonised/GCST90691754.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90691001-GCST90692000/GCST90691754/harmonised/GCST90691754.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90691001-GCST90692000/GCST90691754/harmonised/GCST90691754.h.tsv.gz-meta.yaml


# negative control:
# height (mean, inv-normal transformed)
# 2024 Verma Science study (424,305 European ancestry individuals)
trait=height
mkdir -p ${dir}/${trait}; cd ${dir}/${trait}
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475362/harmonised/GCST90475362.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475362/harmonised/GCST90475362.h.tsv.gz.tbi
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475362/harmonised/GCST90475362.h.tsv.gz-meta.yaml

