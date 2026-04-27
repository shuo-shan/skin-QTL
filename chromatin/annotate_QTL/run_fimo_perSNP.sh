#!/bin/bash
set -euo pipefail

module load bedtools

workingdir=$1
this_snp=$2
prefix=$3

cd "${workingdir}"

genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
genome_fa=/share/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa
memeF=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv12_core_HUMAN_mono_meme_format.meme
export PATH=/home/shuo.shan-umw/meme/bin:/home/shuo.shan-umw/meme/libexec/meme-5.5.5:$PATH

### find SNP row index
rowidx=$(grep -w -n "${this_snp}" "${workingdir}/QTL_${prefix}.bed" | cut -d':' -f1)

if [[ -z "${rowidx}" ]]; then
    echo "ERROR: SNP ${this_snp} not found in ${workingdir}/QTL_${prefix}.bed" >&2
    exit 1
fi

tempdir="${workingdir}/temp${rowidx}"
rm -rf "${tempdir}"
mkdir -p "${tempdir}"
cd "${tempdir}"

### fetch SNP BED row
grep -w "${this_snp}" "${workingdir}/QTL_${prefix}.bed" > SNP.bed

chr=$(awk '{print $1}' SNP.bed)
start=$(awk '{print $2}' SNP.bed)
end=$(awk '{print $3}' SNP.bed)
snp_id=$(awk '{print $4}' SNP.bed)
ref_allele=$(awk '{print $5}' SNP.bed)
alt_allele=$(awk '{print $6}' SNP.bed)

echo "Processing ${snp_id} ${chr}:${start}-${end} REF=${ref_allele} ALT=${alt_allele}"
date

### make 100bp flanking window around SNP: ±50 bp
bedtools slop -b 50 -i SNP.bed -g "${genome}" > SNP_slop100.bed

### get fasta
bedtools getfasta -fi "${genome_fa}" -bed SNP_slop100.bed -fo SNP_slop100.fa -name
echo "got fasta"
date

### build REF and ALT fasta explicitly and compute SNP position within sequence
python3 <<'PY'
# read SNP bed
with open("SNP.bed") as f:
    chrom, start, end, snp_id, ref, alt = f.readline().strip().split()[:6]
start = int(start)
end = int(end)
ref = ref.upper()
alt = alt.upper()

# read slopped window
with open("SNP_slop100.bed") as f:
    wchrom, wstart, wend = f.readline().strip().split()[:3]
wstart = int(wstart)
wend = int(wend)

# SNP position in sequence
snp_idx0 = start - wstart
snp_pos1 = snp_idx0 + 1

# read fasta
lines = [x.strip() for x in open("SNP_slop100.fa") if x.strip()]
header = lines[0]
seq = "".join(lines[1:]).upper()

if snp_idx0 < 0 or snp_idx0 >= len(seq):
    raise ValueError(f"SNP index out of range: snp_idx0={snp_idx0}, seq_len={len(seq)}")

genome_base = seq[snp_idx0]

print(f"[DEBUG] SNP={snp_id}")
print(f"[DEBUG] window={wchrom}:{wstart}-{wend}")
print(f"[DEBUG] SNP genomic={chrom}:{start}-{end}")
print(f"[DEBUG] seq_len={len(seq)}")
print(f"[DEBUG] snp_idx0={snp_idx0}")
print(f"[DEBUG] snp_pos1={snp_pos1}")
print(f"[DEBUG] fasta_base_at_snp={genome_base}")
print(f"[DEBUG] bed_REF={ref}")
print(f"[DEBUG] bed_ALT={alt}")

seq_ref = list(seq)
seq_alt = list(seq)

# force replacement using BED alleles
seq_ref[snp_idx0] = ref
seq_alt[snp_idx0] = alt

with open("SNP_slop100_REF.fa", "w") as out:
    out.write(f">{snp_id}|REF|snp_pos={snp_pos1}\n")
    out.write("".join(seq_ref) + "\n")

with open("SNP_slop100_ALT.fa", "w") as out:
    out.write(f">{snp_id}|ALT|snp_pos={snp_pos1}\n")
    out.write("".join(seq_alt) + "\n")

with open("snp_position_in_seq.txt", "w") as out:
    out.write(str(snp_pos1) + "\n")
PY

snp_pos_seq=$(cat snp_position_in_seq.txt)
echo "SNP position in sequence (1-based) = ${snp_pos_seq}"

### run FIMO on REF and ALT
for allele in REF ALT; do
    fimo --oc "results_${allele}" --thresh 0.001 --verbosity 1 "${memeF}" "SNP_slop100_${allele}.fa"

    ### keep all motif hits except comments/header
    awk 'BEGIN{FS=OFS="\t"} !/^#/ && NR>1' "results_${allele}/fimo.tsv" > "fimo_all_${allele}.tsv"

    ### keep only motifs overlapping the SNP in sequence coordinates
    ### FIMO columns:
    ### 1 motif_id
    ### 2 motif_alt_id
    ### 3 sequence_name
    ### 4 start
    ### 5 stop
    ### 6 strand
    ### 7 score
    ### 8 p-value
    ### 9 q-value
    ### 10 matched_sequence
    awk -v snp="${this_snp}" -v allele="${allele}" -v pos="${snp_pos_seq}" '
        BEGIN{FS=OFS="\t"}
        {
            if ($4 <= pos && $5 >= pos) {
                print snp, allele, $0
            }
        }
    ' "fimo_all_${allele}.tsv" > "fimo_snp_overlap_${allele}.tsv"

done

echo "ran FIMO on REF and ALT"
date

### merge REF/ALT overlapping hits and classify differences
conda activate fastQTL
python3 <<'PY'
import pandas as pd
import numpy as np

def read_or_empty(fn):
    cols = [
        "snp","allele","motif_id","motif_alt_id","sequence_name",
        "start","stop","strand","score","pvalue","qvalue","matched_sequence"
    ]
    try:
        df = pd.read_csv(fn, sep="\t", header=None, names=cols)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=cols)
    return df

ref = read_or_empty("fimo_snp_overlap_REF.tsv")
alt = read_or_empty("fimo_snp_overlap_ALT.tsv")

# use motif_id as primary comparison unit for SNP-overlapping hits
# if multiple hits per motif, keep best p-value
def best_per_motif(df, suffix):
    if df.empty:
        return pd.DataFrame(columns=[
            "motif_id",
            f"start_{suffix}", f"stop_{suffix}", f"strand_{suffix}",
            f"score_{suffix}", f"pvalue_{suffix}", f"qvalue_{suffix}",
            f"matched_sequence_{suffix}"
        ])
    df = df.sort_values(["motif_id", "pvalue"], ascending=[True, True])
    df = df.groupby("motif_id", as_index=False).first()
    return df.rename(columns={
        "start": f"start_{suffix}",
        "stop": f"stop_{suffix}",
        "strand": f"strand_{suffix}",
        "score": f"score_{suffix}",
        "pvalue": f"pvalue_{suffix}",
        "qvalue": f"qvalue_{suffix}",
        "matched_sequence": f"matched_sequence_{suffix}"
    })[[
        "motif_id",
        f"start_{suffix}", f"stop_{suffix}", f"strand_{suffix}",
        f"score_{suffix}", f"pvalue_{suffix}", f"qvalue_{suffix}",
        f"matched_sequence_{suffix}"
    ]]

refb = best_per_motif(ref, "REF")
altb = best_per_motif(alt, "ALT")

merged = refb.merge(altb, on="motif_id", how="outer")

merged["present_REF"] = ~merged["score_REF"].isna()
merged["present_ALT"] = ~merged["score_ALT"].isna()

def classify(row):
    if row["present_REF"] and not row["present_ALT"]:
        return "lost_in_ALT"
    elif not row["present_REF"] and row["present_ALT"]:
        return "gained_in_ALT"
    elif row["present_REF"] and row["present_ALT"]:
        if row["score_ALT"] > row["score_REF"]:
            return "stronger_in_ALT"
        elif row["score_ALT"] < row["score_REF"]:
            return "weaker_in_ALT"
        else:
            return "unchanged"
    else:
        return "absent_both"

merged["delta_score_ALT_minus_REF"] = merged["score_ALT"] - merged["score_REF"]
merged["delta_log10p_ALT_minus_REF"] = (
    -np.log10(merged["pvalue_ALT"]) - (-np.log10(merged["pvalue_REF"]))
)
merged["motif_change_class"] = merged.apply(classify, axis=1)

# optional prioritization score
merged["priority_abs_delta_log10p"] = merged["delta_log10p_ALT_minus_REF"].abs()

merged = merged.sort_values(
    ["motif_change_class", "priority_abs_delta_log10p"],
    ascending=[True, False],
    na_position="last"
)

merged.to_csv("fimo_REF_ALT_comparison.tsv", sep="\t", index=False)
PY

### move final output
mv fimo_REF_ALT_comparison.tsv "${workingdir}/fimo_output_${this_snp}.txt"

cd "${workingdir}"
echo "merged FIMO output written to ${workingdir}/fimo_output_${this_snp}.txt"
