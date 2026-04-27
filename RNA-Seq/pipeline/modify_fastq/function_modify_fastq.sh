#!/bin/bash

fdir=$1
fname=$2
echo "starting to work on" ${fdir}/${fname}; date;
mv ${fdir}/${fname}.p1.fq.gz ${fdir}/archive_${fname}.p1.fq.gz; 
echo "renamed p1..."; date;
zcat ${fdir}/archive_${fname}.p1.fq.gz | awk '{if (NR%4==1) gsub("_",":",$1); print}' | gzip > ${fdir}/${fname}.p1.fq.gz; 
echo "modified p1..."; date;
mv ${fdir}/${fname}.p2.fq.gz ${fdir}/archive_${fname}.p2.fq.gz; 
echo "renamed p2..."; date;
zcat ${fdir}/archive_${fname}.p1.fq.gz | awk '{if (NR%4==1) gsub("_",":",$1); print}' | gzip > ${fdir}/${fname}.p2.fq.gz; 
echo "modified p2..."; date;
rm ${fdir}/archive_${fname}.p1.fq.gz ${fdir}/archive_${fname}.p2.fq.gz; 
echo "done changing seqID for " ${fdir}/${fname}; date;

