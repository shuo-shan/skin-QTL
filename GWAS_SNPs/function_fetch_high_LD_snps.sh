#snp=rsxxxxxx
snp=$1

module load plink
vcfFile=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
#ln -s /pi/manuel.garber-umw/human/skin/eQTLs/fastQTL/fastQTL/GWAS_SNP_with_LD/RegulatoryElementSNPs/plink.37donors.map
#ln -s /pi/manuel.garber-umw/human/skin/eQTLs/fastQTL/fastQTL/GWAS_SNP_with_LD/RegulatoryElementSNPs/plink.37donors.ped
#echo "working on " $snp; date;
plink --file plink.37donors --r2 --ld-snp ${snp} --ld-window-kb 1000 --ld-window-r2 0.8 --out high_ld_variants_of_${snp}
cat high_ld_variants_of_${snp}.ld | awk '{OFS="\t"}NR>1{print "chr"$4,$5,$5+1,$6,$7,$3}' > ${snp}_and_highLD_tags.bed
#cat high_ld_variants_of_${snp}.ld| awk '{OFS="\t"}NR>1{if ($3 != $6) print "chr"$4,$5,$5+1,$6,$7}'
rm high_ld_variants_of_${snp}*; date;

