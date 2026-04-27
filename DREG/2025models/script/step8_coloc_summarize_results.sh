#!/bin/bash
# call Rscript to run coloc on a gene of interest for each condition & QTLtype
# summarize coloc results

# reduces heavy I/O thundering-herd effect in parallel jobs
sleep $((RANDOM % 20))

# ---------- fetch input arguments ------------------- #
ct=$1
trait=$2
cond=$3 # PBS IFNG IFNB TNF
QTLtype=$4 # eQTL or reQTL

DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
echo "starting to perform coloc on all GWAS traits for ${g}"; date


# after jobs run, summarize
mkdir -p ${DIR}/coloc/summary
echo "summarizing ${ct} ${trait} ${cond} ${QTLtype}"; date
cd ${DIR}/coloc/${trait}/${cond}/${QTLtype}
outF=${DIR}/coloc/summary/coloc_${trait}_${cond}_${QTLtype}.txt
rm -f ${outF}
this_f=$(ls coloc_summary_*tsv | head -1)
head -1 ${this_f} > ${outF}
for f in coloc_summary_*.tsv; do
        awk 'NR>1' ${f}  >> ${outF}
done
