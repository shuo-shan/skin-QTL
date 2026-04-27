#!/bin/bash
# make a PDF for gene snp pair of interest given the SNP is in the gene's cis-window.

ct=KRT
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}

# ---------- fetch input arguments ------------------- #
g=$1   #g=ITGA1
snp=$2 #snp=rs2548496
letter=${g:0:1}
mkdir -p ${DIR}/plots/${letter}/${g}
mkdir -p ${DIR}/plots/${letter}/${g}/${snp}
cd ${DIR}/plots/${letter}/${g}/${snp}
echo "starting the script for ${g}:${snp}"; date

# ---------- get gene, SNP, metadata, and modeling stats info ----------- #
# for quick lookup, fetch chunk ID of gene
chunk_dict=${DIR}/data/gene_chunk_dict.txt
chunk_id=$(awk -v gene=${g} '{if ($2==gene) print $1}' ${chunk_dict})

# get gene CPM
cpm_all_f=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt
head -1 ${cpm_all_f} > ${DIR}/plots/${letter}/${g}/${snp}/cpm.txt
awk -v gene=${g} '{OFS=FS="\t"}{if ($NF==gene) print $0}' ${cpm_all_f} >> ${DIR}/plots/${letter}/${g}/${snp}/cpm.txt
echo "got gene CPM"

# get SNP genotype
genotype_chunk_f=${DIR}/chunks/genotype_pairs_chunk_${chunk_id}.tsv
head -1 ${genotype_chunk_f} > ${DIR}/plots/${letter}/${g}/${snp}/genotype.txt
awk -v snp=${snp} '{OFS=FS="\t"}{if ($3==snp) print $0}' ${genotype_chunk_f} >> ${DIR}/plots/${letter}/${g}/${snp}/genotype.txt 
echo "got SNP genotype"

# get model stats
stats_chunk_f=${DIR}/results/result_${chunk_id}.tsv
stats_out_f=${DIR}/plots/${letter}/${g}/${snp}/modelstats.txt
echo -e "celltype\t$(head -1 ${stats_chunk_f})" > ${stats_out_f}
awk -v snp="$snp" -v gene="$g" -v ct="$ct" '
    BEGIN{FS=OFS="\t"}
    NR>1 && $1==snp && $2==gene {print ct, $0}
  ' ${stats_chunk_f} >> ${stats_out_f}
echo "got modeling stats"

# ------------ call Rscript to plot --------------- #
echo "plotting now..."
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${DIR}/step6_make_plot_for_pair_in_window.R ${ct} ${g} ${snp}

cd ${DIR}/plots
rm -r ${DIR}/plots/${letter}/${g}/${snp}
rm -r ${DIR}/plots/${letter}/${g}

echo "done! made plot in ${DIR}/plots/temp_output/plot_${g}_${snp}.pdf"; date










