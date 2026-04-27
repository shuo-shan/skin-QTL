#!/bin/bash
# after running coloc and joining by QTL and GWAS stats, 
# inspect most interesting top genes

# set-up
ct=MEL
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
traits_file=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/traits.txt
gene_file=${DIR}/data/fdr05_genelist.txt


cd ${DIR}

# across all traits, get genes with these characteristics:
# PP.H4 > 0.7 AND p_GWAS < 5e-8 AND p_QTL < 1e-5
# col12           col48             col41
cd ${DIR}/coloc/summary
rm -f top_entries.txt
head -1 coloc_sunburn_IFNB_eQTL.txt > top_entries.txt
for f in coloc_*.txt; do
	awk 'NR>1{if ($12>0.7 && $41<0.00001 && $48<0.00000005) print $0}' ${f} >> top_entries.txt
done
# rank by p_QTL
head -n 1 top_entries.txt > header.txt
tail -n +2 top_entries.txt | sort -k41,41g > body_sorted.txt
cat header.txt body_sorted.txt > top_entries_sorted.txt
mv top_entries_sorted.txt top_entries_PPH4.txt
rm header.txt body_sorted.txt top_entries.txt
res_PPH4=${DIR}/coloc/summary/top_entries_PPH4.txt

# examine main info
cut -f2,3,4,5,12,22,23,27,31,33,50,54,56 ${res_PPH4} | column -t

#
XRRA1	rs11236244	skin_pigmentation
TBL2	rs13233571	skin_pigmentation
UBE2L3	rs3747093	psoriasis
CTSS	rs3768018	skin_pigmentation
CTSS    rs3768018       melanoma dx/hx
MICB	rs6931332	atopic_dermatitis
IRF3	rs8109314	basal_cell_carcinoma
RAC2	rs6000632	vitiligo (Verma 2024)
ZC3H12C	rs199887780	psoriasis
SMARCE1	rs7221109	atopic_dermatitis
SNX32	rs111953392	basal_cell_carcinoma

g=XRRA1
snp=rs11236244
trait=skin_pigmentation
mkdir -p ${DIR}/coloc/summary/best
mkdir -p ${DIR}/coloc/summary/best/${g}
cd ${DIR}/coloc/summary/best/${g}

# [1] Dot plot
bash ${DIR}/step6_plotGeneSnpPair.sh ${g} ${snp}
mv ${DIR}/plots/temp_output/plot_${g}_${snp}.pdf ./

# [2] Coloc manhattan plots
cp ${DIR}/coloc/${trait}/plots/KRT_${g}.locus_tracks.pdf ./
cp ${DIR}/coloc/${trait}/plots_table/KRT_${g}.table.pdf ./

# [3] SNP annotations
bigtable=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_${ct}.txt
head -1 ${bigtable} | tr '\t' '\n' | awk '{print NR,$0}' > rownames.txt
grep ${g} ${bigtable} | grep ${snp} | tr '\t' '\n' > body.txt
if [ -f body.txt ]; then
	echo "yaay!"
	paste rownames.txt body.txt > annotation.txt
	rm body.txt
fi
rm rownames.txt

# [4] Gene annotations
head -1 ${bigtable} > rownames.txt
grep -w ${g} ${bigtable} > body.txt
if [ -f body.txt ]; then
	echo "yaay!"
	cat rownames.txt body.txt > annotation_big.txt
	rm body.txt
fi
rm rownames.txt













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

