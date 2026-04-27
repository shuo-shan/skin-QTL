echo -e "SNP\tCHR\tBP\tP" > all_GWAS123cmh.txt 
for file in GWAS*;do
	echo $file
	awk '{OFS="\t"}NR>1{print $2,$1,$3,$8}' ${file} >> all_GWAS123cmh.txt 
done

