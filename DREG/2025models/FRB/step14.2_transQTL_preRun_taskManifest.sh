#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=1000]"
#BSUB -W 1:00
#BSUB -w "done(transQTL.FRB)"
#BSUB -q short
#BSUB -J makeManifest
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_makeManifest_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/logs/transqtl_makeManifest_%J_%I.err"


ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
taskfile=${DIR}/transQTL/task_manifest.tsv

rm -f "${taskfile}"
touch "${taskfile}"

for cond in PBS IFNB IFNG TNF; do
    for QTLtype in eQTL reQTL; do
        
        # skip impossible combo
        if [[ "${cond}" == "PBS" && "${QTLtype}" == "reQTL" ]]; then
            continue
        fi
        
        chunk_dir=${DIR}/transQTL/chunks/${cond}/${QTLtype}
        
        # skip if folder missing
        [[ -d "${chunk_dir}" ]] || continue
        
        for f in "${chunk_dir}"/gene_QTL_pairs_chunk_*.tsv; do
            [[ -e "$f" ]] || continue
            
            chunk_base=$(basename "$f" .tsv)
            chunk_id=$(echo "$chunk_base" | sed 's/gene_QTL_pairs_chunk_//')
            
            printf "%s\t%s\t%s\t%s\n" "${cond}" "${QTLtype}" "${chunk_id}" "${f}" >> "${taskfile}"
        done
    done
done

echo "Task file written to: ${taskfile}"
wc -l "${taskfile}"
