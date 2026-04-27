#############################################################
cd /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20211025_Celseq2
f=files_for_function_run_ESAT.txt
outdir=$(head -1 ${f} | cut -d$'\t' -f3)
cd ${outdir}
mkdir scripture
mkdir gene
mkdir window
mkdir umi_distributions
cd /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20211025_Celseq2
while read line;do
  param1=$(echo ${line} | cut -d' ' -f1)
  param2=$(echo ${line} | cut -d' ' -f2)
  param3=$(echo ${line} | cut -d' ' -f3)
  prefix=$(echo ${param2} | sed "s/.bam//g")
  echo "bash function_run_ESAT.sh" ${param1} ${param2} ${param3} | bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=110000] -R "select[rh=6]" -q long -J ${prefix}  -o "${param3}/../log/%J%I.ESAT.${prefix}.out" -e "${param3}/../log/%J%I.ESAT.${prefix}.err"
done < ${f}

##############################################################
## organize into one
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/ESAT/output/downsample/gene
#cd $dir
## paste is appropriate here b/c all files have the same row names
#paste * > temp.gene.txt
## get rid of duplicated gene symbol, chr, strand columns from pasting:
#keep=$(head -1 temp.gene.txt | tr '\t' '\n' | grep -v -n -E 'Symbol|chr|strand' | cut -d":" -f1 | tr '\n' ',' | sed 's/\(.*\),/\1/' | sed 's/4,/1,2,3,4,/')
## get rid of sample barcode from sample name
#cat temp.gene.txt | cut -f$keep | sed '1 s/:[ATCG]\{6\}//g' > all.gene.txt
#rm temp.gene.txt
#chmod 777 all.gene.txt
