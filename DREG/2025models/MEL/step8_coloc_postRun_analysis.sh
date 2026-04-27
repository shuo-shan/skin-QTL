#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 08:00
#BSUB -J plot_coloc_gene_MEL
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/plot_coloc_genes_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/plot_coloc_genes_%J.err"
# after running coloc and joining by QTL and GWAS stats, 
# inspect most interesting top genes

# set-up
ct=MEL
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
traits_file=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/traits.txt
gene_file=${DIR}/data/fdr05_genelist.txt
cd ${DIR}


# go to Rscript to run: step8_coloc_postRun_analysis_annotation_explorer.R
# this outputs files in coloc/summary/master_coloc_trait_summary.tsv
# sig gene:snp pairs are in that file

# get all gene SNP pairs that coloc (PPH4>0.7)
f=${DIR}/coloc/summary/master_coloc_trait_summary.tsv
while read trait;do
  awk 'NR>1' ${f} | grep ${trait} | cut -f10 | tr ',' '\n' | grep ":" > ${DIR}/coloc/summary/${trait}_eGenes_leadQTL_pairs.txt
  echo ${trait}
  while read line;do
	  g=$(echo ${line} | cut -d':' -f1)
	  snp=$(echo ${line} | cut -d':' -f2)

	  echo ${f} ${snp}
	  mkdir -p ${DIR}/coloc/summary/best
	  mkdir -p ${DIR}/coloc/summary/best/${g}
	  cd ${DIR}/coloc/summary/best/${g}

	  # store which trait & QTL condition x QTLtype is this gene colocalizing in
	  awk 'NR>1' ${f} | grep ${trait} | grep "${g}:" | cut -f1,2,3,4 >> coloc_annotations.txt

	  # [1] Dot plot
	  bash ${DIR}/step6_plotGeneSnpPair.sh ${g} ${snp}
	  mv ${DIR}/plots/temp_output/plot_${g}_${snp}.pdf ./

	  # [2] Coloc manhattan plots
	  cp ${DIR}/coloc/${trait}/plots/${ct}_${g}.locus_tracks.pdf ./
	  cp ${DIR}/coloc/${trait}/plots_table/${ct}_${g}.table.pdf ./

	  # [3] SNP annotations
	  bigtable=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_${ct}.txt
          head -1 ${bigtable} | tr '\t' '\n' | awk '{print NR,$0}' > rownames.txt
          grep ${g} ${bigtable} | grep ${snp} | tr '\t' '\n' > body.txt
          if [ -f body.txt ]; then
                  echo "yaay!"
                  paste rownames.txt body.txt > annotation_${snp}.txt
                  rm body.txt
          fi
          rm rownames.txt
          
          # [4] Gene annotations
          head -1 ${bigtable} > rownames.txt
          grep -w ${g} ${bigtable} > body.txt
          if [ -f body.txt ]; then
                  echo "yaay!"
                  cat rownames.txt body.txt > annotation_big_${g}.txt
                  rm body.txt
          fi
          rm rownames.txt
  done < ${DIR}/coloc/summary/${trait}_eGenes_leadQTL_pairs.txt          
done < ${traits_file}



## example
#g=XRRA1
#snp=rs11236244
#trait=skin_pigmentation
#mkdir -p ${DIR}/coloc/summary/best
#mkdir -p ${DIR}/coloc/summary/best/${g}
#cd ${DIR}/coloc/summary/best/${g}
#
## [1] Dot plot
#bash ${DIR}/step6_plotGeneSnpPair.sh ${g} ${snp}
#mv ${DIR}/plots/temp_output/plot_${g}_${snp}.pdf ./
#
## [2] Coloc manhattan plots
#cp ${DIR}/coloc/${trait}/plots/MEL_${g}.locus_tracks.pdf ./
#cp ${DIR}/coloc/${trait}/plots_table/MEL_${g}.table.pdf ./
#
## [3] SNP annotations
#bigtable=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_${ct}.txt
#head -1 ${bigtable} | tr '\t' '\n' | awk '{print NR,$0}' > rownames.txt
#grep ${g} ${bigtable} | grep ${snp} | tr '\t' '\n' > body.txt
#if [ -f body.txt ]; then
#	echo "yaay!"
#	paste rownames.txt body.txt > annotation.txt
#	rm body.txt
#fi
#rm rownames.txt
#
## [4] Gene annotations
#head -1 ${bigtable} > rownames.txt
#grep -w ${g} ${bigtable} > body.txt
#if [ -f body.txt ]; then
#	echo "yaay!"
#	cat rownames.txt body.txt > annotation_big.txt
#	rm body.txt
#fi
#rm rownames.txt













# ---------------------------------------------------------
# across all traits, get genes with these characteristics:
# PP.H3 > 0.8 AND p_GWAS < 5e-8 AND p_QTL < 1e-6
# col11           col48             col41
cd ${DIR}/coloc/summary
rm -f top_entries.txt
head -1 coloc_sunburn_IFNB_eQTL.txt > top_entries.txt
for f in coloc_*.txt; do
        awk 'NR>1{if ($11>0.8 && $41<0.000001 && $48<0.00000005) print $0}' ${f} >> top_entries.txt
done
# rank by p_QTL
head -n 1 top_entries.txt > header.txt
tail -n +2 top_entries.txt | sort -k41,41g > body_sorted.txt
cat header.txt body_sorted.txt > top_entries_sorted.txt
mv top_entries_sorted.txt top_entries_PPH3.txt
rm header.txt body_sorted.txt top_entries.txt

# examine main info
cut -f2,3,4,5,12,22,23,27,31,33,50,54,56 ${DIR}/coloc/summary/top_entries_PPH3.txt | column -t

