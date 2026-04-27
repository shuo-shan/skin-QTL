#!/bin/bash
# summarize the number of genes with a eQTL or reQTL in each condition

DIR="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/"
OUTDIR=${DIR}/data

# ----------------------------------------------------------------- #
# tabularize gene names
OUTF=${OUTDIR}/all_QTL_genes_FDR05.txt
echo -e "celltype\tcondition\tQTLtype\tgene\tlead_snp\tpmin\tq_gene" > ${OUTF}
cat ${DIR}/MEL/eigenMT/results/*_gene_fdr05_table.txt | awk '{OFS=FS="\t"}{print $1,$2,$3,$4,$8,$7,$11}' | grep -v "lead_snp" >> ${OUTF}
cat ${DIR}/KRT/eigenMT/results/*_gene_fdr05_table.txt | awk '{OFS=FS="\t"}{print $1,$2,$3,$4,$9,$8,$13}' | grep -v "lead_snp" >> ${OUTF}
cat ${DIR}/FRB/eigenMT/results/*_gene_fdr05_table.txt | awk '{OFS=FS="\t"}{print $1,$2,$3,$4,$9,$8,$13}' | grep -v "lead_snp" >> ${OUTF}


# summarize number of discoveries
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/script/step9_summarize_all_QTL_genes.R

