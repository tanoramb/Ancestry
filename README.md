# Ancestry

Steps for performing ancestry analisis using RFMix 1.5 or 2.0 and ShapeIT or Beagle

## 1. Data preparation

We have to process original data (VCF files) by removing duplications, and we do it for each chromosome in a parallel fashion.

Links:

+ [TABIX 0.2.6](https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2/download)
+ [BCFTOOLS 1.10.2](https://sourceforge.net/projects/samtools/files/samtools/1.10.2/bcftools-1.10.2.tar.bz2/download)
+ [SCRIPTS](https://github.com/tanoramb/Ancestry/tree/eaf4b0bf5a445748de9eccc83b78c008d61b2b34/scripts)



```
#! /bin/bash
RMD="/path/to/script/remove_duplicate_by_pos.py"
TABIX="/path/to/software/Tabix/tabix-0.2.6/tabix"
PARALLEL="/path/to/script/parallel.sh"

NCPU=10 #Number of CPUs

# VCF files are named as "1KG_PEL_chr*.recode.vcf" where "*" is the number of the chromosome 1 through 22
# Then the VCF files are processed to remove duplicated positions, bgzipped, and indexed.
${PARALLEL} -j ${NCPU} -r "${RMD} -i 1KG_PEL_chr*.recode.vcf -o 1KG_PEL_chr*.recode.rmd.vcf > rmd.chr*.log 2>&1" $(seq 1 22)
${PARALLEL} -j ${NCPU} -r "bgzip 1KG_PEL_chr*.recode.rmd.vcf" $(seq 1 22)
${PARALLEL} -j ${NCPU} -r "${TABIX} -p vcf 1KG_PEL_chr*.recode.rmd.vcf.gz" $(seq 1 22)
```

If the VCF files must be further processed to remove/keep samples:

```
BCFTOOLS="/path/to/software/BCFTools/bcftools-1.10.2/bcftools"

# VCF files are gunzipped
${PARALLEL} -j ${NCPU} -r "gunzip -k 1KG_PEL_chr*.recode.rmd.vcf.gz" $(seq 1 22)

IDS90="ids_pel_ge90.txt" #Text file containing Sample IDs to keep from VCF file. In this case, those Samples with more than 90% Native-American ancestry (analysis done elsewhere)
IDS95="ids_pel_ge95.txt" #Text file containing Sample IDs to keep from VCF file. In this case, those Samples with more than 95% Native-American ancestry (analysis done elsewhere)

${PARALLEL} -j ${NCPU} -r "gunzip -k 1KG_PEL_chr*.recode.rmd.vcf.gz" $(seq 1 22)
for PCT in 90 95; do
    IDS="ids_pel_ge${PCT}.txt"
    ${PARALLEL} -j ${NCPU} -r "${BCFTOOLS} view -S ${IDS} 1KG_PEL_chr*.recode.rmd.vcf > 1KG_PEL_chr*.recode.rmd.ge${PCT}.vcf" $(seq 1 22)
    ${PARALLEL} -j ${NCPU} -r "bgzip 1KG_PEL_chr*.recode.rmd.ge${PCT}.vcf" $(seq 1 22)
    ${PARALLEL} -j ${NCPU} -r "${TABIX} -p vcf 1KG_PEL_chr*.recode.rmd.ge${PCT}.vcf.gz" $(seq 1 22)
done
\rm -rf .tmp* *.vcf

```
