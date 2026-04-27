

conda activate HLAtyping



###### 08192024
dir=/pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing
cd ${dir}/data
### retrieve the chr6 of the genotyping file.
module load bcftools
bcftools index bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.vcf.gz 
bcftools view -r chr6 bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.vcf.gz -Oz -o bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.chr6.vcf.gz
vcfF=/pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/data/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.chr6.vcf.gz

module load plink
plink --vcf ${vcfF} --make-bed --out skineQTL_chr6

######## QC of genotype data ######## 
# https://github.com/immunogenomics/HLA_analyses_tutorial/blob/main/tutorial_HLAQCImputation.ipynb
# overview: 3,752,078 variants on chr6, 37 donors.

#### 1. dedup and remove high missing rate variants.
# plink flags duplicated variants and high missing rate by position (via -missing command)
plink --bfile skineQTL_chr6 --missing --out skineQTL_chr6
# SNP2HLA python script excludes those flagged variants.
python ${dir}/scripts/get_duprem_var.py skineQTL_chr6
plink --bfile skineQTL_chr6 --exclude skineQTL_chr6.remdup.snp --make-bed --out skineQTL_chr6.dedup
#### conclusion: removed 0 SNP.


#### reverse/forward strand flip
module load gcc/12.2.0
module load python2
pip install snpflip
conda activate snpflip
snpflip --fasta-genome=/share/data/umw_biocore/genome_data/human/hg38/hg38.fa --bim-file=skineQTL_chr6.dedup.bim --output-prefix=skineQTL_chr6.dedup.test
plink --bfile skineQTL_chr6.dedup --exclude skineQTL_chr6.dedup.test.ambiguous --make-bed --out skineQTL_chr6.dedup.ambstrandrem
module unload gcc/12.2.0 python2
#### conclusion: removed 0 SNP.

#### Remove palindromic SNPs
# Remove A/T C/G SNPs
awk '(($5=="C"&&$6=="G")||($5=="G"&&$6=="C")||($5=="A"&&$6=="T")||($5=="T"&&$6=="A")) {print $2} ' skineQTL_chr6.dedup.ambstrandrem.bim > AT_CG_SNPS.txt
plink --bfile skineQTL_chr6.dedup.ambstrandrem --exclude AT_CG_SNPS.txt --make-bed --out skineQTL_chr6.dedup.ambstrandrem
#### conclusion: removed 0 SNP.

#### Remove poor quality variants (variant missningness of 10% and hwe 1e-10)
plink --bfile skineQTL_chr6.dedup.ambstrandrem --geno 0.1 --hwe 1e-10 --make-bed --out skineQTL_chr6.dedup.ambstrandrem.1stSNPQC
#### conclusion: 67895 variants removed due to missing genotype data (--geno). --hwe: 0 variants removed due to Hardy-Weinberg exact test.

#### skipped: sample-level QC and other QC steps.


#### 2. extract the MHC region: final input for imputation
#isolate MHC
plink --bfile skineQTL_chr6.dedup.ambstrandrem.1stSNPQC --chr 6 --from-mb 28 --to-mb 34 --make-bed --out skineQTL_chr6.final
# 106455 out of 3091473 variants loaded from .bim file

#### 3. rename variants for SNP2HLA input format
mv skineQTL_chr6.final.bim skineQTL_chr6.final.bim.old
python ${dir}/scripts/rename_bim.py ${dir}/data/Tutorial_1KGonly.bim skineQTL_chr6.final.bim

#### 4. Genotype phasing required for HLA imputaion
shapeitDir=/pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/tutorial/HLA_analyses_tutorial/SHAPEIT/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin
${shapeitDir}/shapeit -B skineQTL_chr6.final -M ${dir}/data/genetic_map_chr6_combined_b37.txt -O skineQTL_chr6.final.shapeit.phased --thread 8 --seed 0 --output-log skineQTL_chr6.final.shapeit.phased.log
${shapeitDir}/shapeit -convert --input-haps skineQTL_chr6.final.shapeit.phased --output-vcf skineQTL_chr6.final.shapeit.phased.vcf
conda activate fastQTL
bgzip -c skineQTL_chr6.final.shapeit.phased.vcf > skineQTL_chr6.final.shapeit.phased.vcf.gz
tabix skineQTL_chr6.final.shapeit.phased.vcf.gz
conda deactivate

#### 5. HLA imputation
plink --vcf skineQTL_chr6.final.shapeit.phased.vcf.gz --keep-allele-order --make-bed --out skineQTL_chr6.final.shapeit.phased
plink --bfile skineQTL_chr6.final.shapeit.phased --freq --out skineQTL_chr6.final.shapeit.phased.FRQ
zcat skineQTL_chr6.final.shapeit.phased.vcf.gz | grep -v "#" | awk '{print $3,$2,$4,$5}' > Tutorial_1KGonly.bgl.phased.markers
snp2hlaDir=/pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/tutorial/HLA_analyses_tutorial
cd ${snp2hlaDir}/scripts

python ./SNP2HLA.py \
  -i hgdp_chr6.final \
  -o hgdp_chr6.final.SNP2HLApy.imputed \
  -rf Tutorial_1KGonly \
  --nthreads 10 \
  --java-mem=80g --tolerated-diff=0.5

python ./SNP2HLA.py \
	-i ${dir}/data/skineQTL_chr6.final \
	-o ${dir}/data/skineQTL_chr6.final.SNP2HLApy.imputed \
	-rf ${dir}/data/skineQTL_chr6.final.shapeit.phased \
	--nthreads 10 \
	--java-mem=80g --tolerated-diff=0.5


java -Djava.io.tmpdir=/pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/output/skineQTL_chr6.final.SNP2HLApy.imputed.javatmpdir \
        -Xmx80000m \
        -jar ./beagle.jar \
        gt=/pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/data/skineQTL_chr6.final.shapeit.phased.vcf.gz \
        impute=true gprobs=true nthreads=10 chrom=6 niterations=5 lowmem=true  \
        out=/pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/output/skineQTL_chr6.final.SNP2HLApy.imputed.bgl.phased



./plink --silent --allow-no-sex --make-bed --vcf /pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/output/skineQTL_chr6.final.SNP2HLApy.imputed.bgl.phased.vcf.gz --a1-allele /pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/data/skineQTL_chr6.final.shapeit.phased.markers 4 1 --out /pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/output/skineQTL_chr6.final.SNP2HLApy.imputed



        # (2) Imputation result in *.{bed,bim,fam} files (*.vcf.gz => *.{bed,bim,fam})
        command = ' '.join([PLINK, "--make-bed", "--vcf", __IMPUTED__, "--a1-allele {} 4 1".format(_reference_panel+".markers"), "--out", OUTPUT])
        #print(command)
        #os.system(command)


        # (3) Dosage file
        command = ' '.join(["gunzip -c", __IMPUTED__, "|", "cat", "|", "java -jar {} > {}".format(_vcf2gprobs, OUTPUT+".bgl.gprobs")])
        #print(command)
        #os.system(command)

        __gprobs__ = OUTPUT+".bgl.gprobs"


        command = ' '.join(["tail -n +2 {}".format(__gprobs__), "|",
                            PARSEDOSAGE, "- > {}".format(OUTPUT+".dosage")])
        #print(command)
        #os.system(command)


cd ${dir}/output
bcftools view -H skineQTL_chr6.final.SNP2HLApy.imputed.bgl.phased.vcf.gz | grep rs | cut -f3 | grep rs | sort | uniq -c | tr -s " " | awk '{if ($1>1) print $2}' > duplicate_variants.txt
mv ${dir}/output/duplicate_variants.txt ${dir}/data/duplicate_variants.txt
cd ${dir}/data
bcftools view --exclude ID=@duplicate_variants.txt bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.chr6.vcf.gz -Oz -o bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.chr6.dedup.vcf.gz
vcfF=bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.chr6.dedup.vcf.gz
module load plink
plink --vcf ${vcfF} --make-bed --out skineQTL_chr6

#### Remove poor quality variants (variant missningness of 10% and hwe 1e-10)
plink --bfile skineQTL_chr6 --geno 0.1 --hwe 1e-10 --make-bed --out skineQTL_chr6.dedup.ambstrandrem.1stSNPQC
#### conclusion: 81,548 variants removed due to missing genotype data (--geno). --hwe: 0 variants removed due to Hardy-Weinberg exact test.
#### 3,669,667 variants and 37 people pass filters and QC.
#### skipped: sample-level QC and other QC steps.

#### 2. extract the MHC region: final input for imputation
#isolate MHC
plink --bfile skineQTL_chr6.dedup.ambstrandrem.1stSNPQC --chr 6 --from-mb 28 --to-mb 34 --make-bed --out skineQTL_chr6.final
# 106455 out of 3091473 variants loaded from .bim file

#### 3. rename variants for SNP2HLA input format
mv skineQTL_chr6.final.bim skineQTL_chr6.final.bim.old
python ${dir}/scripts/rename_bim.py ${dir}/data/Tutorial_1KGonly.bim skineQTL_chr6.final.bim

#### 4. Genotype phasing required for HLA imputaion
shapeitDir=/pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/tutorial/HLA_analyses_tutorial/SHAPEIT/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin
${shapeitDir}/shapeit -B skineQTL_chr6.final -M ${dir}/data/genetic_map_chr6_combined_b37.txt -O skineQTL_chr6.final.shapeit.phased --thread 8 --seed 0 --output-log skineQTL_chr6.final.shapeit.phased.log
${shapeitDir}/shapeit -convert --input-haps skineQTL_chr6.final.shapeit.phased --output-vcf skineQTL_chr6.final.shapeit.phased.vcf
conda activate fastQTL
bgzip -c skineQTL_chr6.final.shapeit.phased.vcf > skineQTL_chr6.final.shapeit.phased.vcf.gz
tabix skineQTL_chr6.final.shapeit.phased.vcf.gz
conda deactivate

#### 5. HLA imputation
plink --vcf skineQTL_chr6.final.shapeit.phased.vcf.gz --keep-allele-order --make-bed --out skineQTL_chr6.final.shapeit.phased
plink --bfile skineQTL_chr6.final.shapeit.phased --freq --out skineQTL_chr6.final.shapeit.phased.FRQ
zcat skineQTL_chr6.final.shapeit.phased.vcf.gz | grep -v "#" | awk '{print $3,$2,$4,$5}' > Tutorial_1KGonly.bgl.phased.markers
snp2hlaDir=/pi/manuel.garber-umw/human/skin/eQTLs/HLA_typing/tutorial/HLA_analyses_tutorial
cd ${snp2hlaDir}/scripts

python ./SNP2HLA.py \
  -i hgdp_chr6.final \
  -o hgdp_chr6.final.SNP2HLApy.imputed \
  -rf Tutorial_1KGonly \
  --nthreads 10 \
  --java-mem=80g --tolerated-diff=0.5

python ./SNP2HLA.py \
        -i ${dir}/data/skineQTL_chr6.final \
        -o ${dir}/data/skineQTL_chr6.final.SNP2HLApy.imputed \
        -rf ${dir}/data/skineQTL_chr6.final.shapeit.phased \
        --nthreads 10 \
        --java-mem=80g --tolerated-diff=0.5


