# Ancestry

Steps for performing ancestry analisis using RFMix 1.5 or 2.0 and ShapeIT or Beagle

## 1. Data preparation

We have to process original data (VCF files) by removing duplications, and we do it for each chromosome in a parallel fashion:

```
#! /bin/bash
RMD="/path/to/script/remove_duplicate_by_pos.py"
TABIX="/path/to/software/Tabix/tabix-0.2.6/tabix"
PARALLEL="/path/to/script/parallel.sh"

NCPU=22

${PARALLEL} -j ${NCPU} -r "${RMD} -i 1KG_PEL_chr*.recode.vcf -o 1KG_PEL_chr*.recode.rmd.vcf > rmd.chr*.log 2>&1" $(seq 1 22)
${PARALLEL} -j ${NCPU} -r "bgzip 1KG_PEL_chr*.recode.rmd.vcf" $(seq 1 22)
${PARALLEL} -j ${NCPU} -r "${TABIX} -p vcf 1KG_PEL_chr*.recode.rmd.vcf.gz" $(seq 1 22)
${PARALLEL} -j ${NCPU} -r "gunzip -k 1KG_PEL_chr*.recode.rmd.vcf.gz" $(seq 1 22)
```
