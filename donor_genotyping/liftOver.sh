java -jar /share/pkg/picard/2.23.3/picard.jar CreateSequenceDictionary R=hg38.fa O=hg38.dict
grep ">" hg38.fa  | sed 's/>//g' | awk '{print $1"\t"$1}' | sed 's/chr//1' > chr_name_conv.txt
bcftools annotate --rename-chrs chr_name_conv.txt 481d8927-d82c-4865-b00b-d530f346041c_impute-vcf-merged.vcf.bgz | bgzip > renamed_impute-vcf-merged.vcf.bgz
tabix -p vcf renamed_impute-vcf-merged.vcf.bgz 
