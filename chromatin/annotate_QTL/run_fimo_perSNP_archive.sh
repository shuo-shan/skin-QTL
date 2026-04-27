#!/bin/bash
module load bedtools

workingdir=$1
this_snp=$2 # SNP ID
prefix=$3

cd ${workingdir}

### go to temp working folder
rowidx=$(grep -w -n ${this_snp} ${workingdir}/QTL_${prefix}.bed | cut -d':' -f1)
rm -r temp${rowidx}
mkdir temp${rowidx}
cd temp${rowidx}

### fetch SNP bed
grep -w ${this_snp} ${workingdir}/QTL_${prefix}.bed > SNP.bed

### Slop to 100bp
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
bedtools slop -b 50 -i SNP.bed -g ${genome} > SNP_slop100.bed

### getfasta (REF)
genome_fa=/share/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa
bedtools getfasta -fi ${genome_fa} -bed SNP_slop100.bed -fo SNP_slop100.fa
echo "got fasta";date

### run FIMO on 100slop region flanking SNP in SNP=REF and SNP=ALT two scenarios
export PATH=/home/shuo.shan-umw/meme/bin:/home/shuo.shan-umw/meme/libexec/meme-5.5.5:$PATH
memeF=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv12_core_HUMAN_mono_meme_format.meme

ref_allele=$(awk '{print $5}' SNP.bed)
alt_allele=$(awk '{print $6}' SNP.bed)
snp_pos=$(awk '{print $2}' SNP.bed)

for allele in REF ALT; do
    if [ "${allele}" == "REF" ]; then
        this_allele=${ref_allele}
    else
        this_allele=${alt_allele}
    fi

    python3 -c "
lines = open('SNP_slop100.fa').readlines()
seq = list(lines[1].strip())
seq[50] = '${this_allele}'
print(lines[0].strip())
print(''.join(seq))
" > SNP_slop100_${allele}.fa

    fimo --oc results_${allele} --thresh 0.001 --verbosity 1 ${memeF} SNP_slop100_${allele}.fa

    # at this step, filter for motif that overlap SNP, then get top100
    output_file="fimo_SNPoverlapping_top100_${allele}.txt"
    rm -f ${output_file}
    cat results_${allele}/fimo.tsv | \
        grep .H12CORE. | \
        awk -v pos=${snp_pos} '$6 <= pos && $7 > pos' | \
        awk -v snp=${this_snp} -v allele=${allele} 'BEGIN{OFS="\t"}{print snp,allele,$0}' \
        >> ${output_file}

    # get all motifs in this 100bp region, top 100
    output_file=fimo_top100_${allele}.txt
    rm -f ${output_file}
    cat results_${allele}/fimo.tsv | \
        grep .H12CORE. | \
        awk -v snp=${this_snp} -v allele=${allele} 'BEGIN{OFS="\t"}{print snp,allele,$0}' \
        >> ${output_file}

done
echo "ran FIMO on SNP's REF and ALT alleles";date

# if I want top 100
#    cat results_${allele}/fimo.tsv | \
#        grep .H12CORE. | \
#        awk -v pos=${snp_pos} '$6 <= pos && $7 > pos' | \
#        head -100 | \
#        awk -v snp=${this_snp} -v allele=${allele} 'BEGIN{OFS="\t"}{print snp,allele,$0}' \
#        >> ${output_file}


### merge tables
# merge REF and ALT top100, annotate motif_overlap_SNP_REF and motif_overlap_SNP_ALT
output_merged="fimo_top100_merged.txt"
rm -f ${output_merged}

# get all unique motifs across REF and ALT
cat fimo_top100_REF.txt fimo_top100_ALT.txt | \
    awk 'BEGIN{OFS="\t"} {
        key=$1"\t"$3"\t"$4
        if (!(key in seen)) {
            seen[key] = 1
            lines[key] = $0
        }
    } END {
        for (key in lines) print lines[key]
    }' | \
    awk -v ref_file="fimo_SNPoverlapping_top100_REF.txt" \
        -v alt_file="fimo_SNPoverlapping_top100_ALT.txt" \
    'BEGIN{
        OFS="\t"
        while((getline line < ref_file) > 0) { split(line,a,"\t"); ref_motifs[a[3]] = 1 }
        while((getline line < alt_file) > 0) { split(line,a,"\t"); alt_motifs[a[3]] = 1 }
	print "SNP\tallele\tmotif_id\tchr\tstart\tend\tstrand\tscore\tpval\tqval\tsequence\tmotif_overlap_SNP_REF\tmotif_overlap_SNP_ALT"
    }
    {
        ref_flag = ($3 in ref_motifs) ? "TRUE" : "FALSE"
        alt_flag = ($3 in alt_motifs) ? "TRUE" : "FALSE"
        print $0, ref_flag, alt_flag
    }' \
    >> ${output_merged}



### write to file
mv ${output_merged} ${workingdir}/fimo_output_${this_snp}.txt
cd ${workingdir}
#rm -r temp${rowidx}
echo "merged fimo output written to ${workingdir}/fimo_output_${this_snp}.txt"

#### acquire lock before writing
#output_file="${workingdir}/fimo_output.txt"
#lock_file="${workingdir}/fimo_output.lock"
#exec 9>"$lock_file" # opens the lock file for writing and associates it with file descriptor 9. lock file is created if it doesn't already exist.
#flock 9 # apply an advisory lock on this file
## write to the output file
#cat ${output_merged} | awk -v snp=${this_snp} '{OFS=FS="\t"}{print snp,$0}' >> ${output_file}
#cat ${output_merged} > ${workingdir}/
## release the lock
#flock -u 9

