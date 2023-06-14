# Ancestry Pipeline

Steps for performing ancestry analisis using RFMix 1.5 or 2.0 and ShapeIT or Beagle

## 1. VCF Pre-processing

We have to process original data (VCF files) by removing duplications, and we do it for each chromosome in a parallel fashion.

Links:

+ [TABIX 0.2.6](https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2/download)
+ [BCFTOOLS 1.10.2](https://sourceforge.net/projects/samtools/files/samtools/1.10.2/bcftools-1.10.2.tar.bz2/download)
+ [SCRIPTS](https://github.com/tanoramb/Ancestry/tree/eaf4b0bf5a445748de9eccc83b78c008d61b2b34/scripts)



```
#! /bin/bash
RMD="/path/to/script/remove_duplicate_by_pos.py" #e.g. /data/tanoramb/Scripts/v2/remove_duplicate_by_pos.py in Astro
TABIX="/path/to/software/Tabix/tabix-0.2.6/tabix"
PARALLEL="/path/to/script/parallel.sh" #e.g. /data/tanoramb/Scripts/other/parallel.sh in Astro 

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

IDS90="ids_pel_ge90.txt" #Text file containing Sample IDs to keep from VCF file. In this case, those Samples with more than 90% Native-American ancestry (analysis done elsewhere)
IDS95="ids_pel_ge95.txt" #Text file containing Sample IDs to keep from VCF file. In this case, those Samples with more than 95% Native-American ancestry (analysis done elsewhere)

# VCF files are gunzipped
${PARALLEL} -j ${NCPU} -r "gunzip -k 1KG_PEL_chr*.recode.rmd.vcf.gz" $(seq 1 22)

# We run bcftools for both percentages 90% and 95%, and to keep the corresponding samples. Then the files are bgzipped and indexed.
for PCT in 90 95; do
    IDS="ids_pel_ge${PCT}.txt"
    ${PARALLEL} -j ${NCPU} -r "${BCFTOOLS} view -S ${IDS} 1KG_PEL_chr*.recode.rmd.vcf > 1KG_PEL_chr*.recode.rmd.ge${PCT}.vcf" $(seq 1 22)
    ${PARALLEL} -j ${NCPU} -r "bgzip 1KG_PEL_chr*.recode.rmd.ge${PCT}.vcf" $(seq 1 22)
    ${PARALLEL} -j ${NCPU} -r "${TABIX} -p vcf 1KG_PEL_chr*.recode.rmd.ge${PCT}.vcf.gz" $(seq 1 22)
done
\rm -rf .tmp* *.vcf

```

## 2. Data preparation

In this step, the data will be prepared to calculate the local ancestry later. The steps described in the master script (per chromosome) involve:

+ Defining variables to store paths to scripts and software.
+ Defining variables to store paths to relevant input data.
+ Getting VCF files containing reference positions (in this case, Chilean positions).
+ Intersecting all VCF files and correcting genotypes (e.g. strand issues).
+ Phasing the intersected/corrected genotypes (using Beagle5 or ShapeIT).
+ Generating proper input files for RFMix 1.5 or RFMix 2.0.

Links:

+ [BEAGLE 5.4](http://faculty.washington.edu/browning/beagle/beagle.html)
+ [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
+ [PLINK 1](https://zzz.bwh.harvard.edu/plink/download.shtml#download)
+ [PLINK 2](https://www.cog-genomics.org/plink/2.0/)
+ [VCFTOOLS](https://vcftools.github.io/man_latest.html)
+ [SCRIPTS](https://github.com/tanoramb/Ancestry/tree/c0557df321d6513d211e4993f941c39fe6e1a57a/scripts)

```
#! /bin/bash
# Script that calls the master script for all chromosomes in a parallel fashion for each chromosome
  
NCPU=10 #Number of CPUs

PARALLEL="/path/to/script/parallel.sh" #e.g. /data/tanoramb/Scripts/other/parallel.sh in Astro
SCRIPT="/path/to/script/master_script_data_beagle5_rfmix15and20.sh"
#SCRIPT="/path/to/script/master_script_data_shapeit_rfmix20.sh" #If the phasing must be done using ShapeIT and RFMix 2.0 only
${PARALLEL} -j ${NCPU} -r "${SCRIPT} -c=* > data.chr*.log 2>&1" 1 22 2 21 3 20 4 19 5 18 6 17 7 16 8 15  9 14 10 13 11 12

#IF RFMix 2.0 is run later
POPS=("PEL" "IBS" "YRI")
echo "" > refs.sample_map.tsv; :> refs.sample_map.tsv
for POP in ${POPS[*]}; do
    name=$(echo ${POP} | awk '{print tolower($0)}')
    gunzip -c phased.intersected.norm.${name}.chr22.vcf.gz | grep "#CHROM" | cut -f10- | sed 's/\t/\n/g' | awk -v group="${POP}" '{print $1"\t"group}' >> refs.sample_map.tsv
done
```


## 3. Ancestry calculation

In this step, RFMix is called to calculate ancestries.

For the case of RFMix 2.0:

```
#! /bin/bash
  
NCPU=10 #Number of CPUs

PARALLEL="/path/to/script/parallel.sh" #e.g. /data/tanoramb/Scripts/other/parallel.sh in Astro

RFMIX20="/path/to/Software/RFmix/RFMix2/rfmix"  #e.g. /data/tanoramb/Software/RFmix/RFMix2/rfmix in Astro
MAPS="/path/to/Data/recombinationmaps/rfmix" #e.g. /data/tanoramb/Data/recombinationmaps/rfmix in Astro
DIR="/path/to/data" #Path to directory containing files and subdirectories generated by Master script of previous step 

${PARALLEL} -j ${NCPU} -r "${RFMIX20} -f ${DIR}/phased.intersected.norm.chi.chr*.vcf.gz -r ${DIR}/phased.intersected.norm.REFS.chr*.vcf.gz -m ${DIR}/refs.sample_map.tsv -g ${MAPS}/chr*.rfmix.map -o rfmix20.chi.chr*.deconvoluted --chromosome=* > rfmix20.chr*.log 2>&1" 1 22 2 21 3 20 4 19 5 18 6 17 7 16 8 15 9 14 10 13 11 12
```

For the case of RFMix 1.5:

```
#! /bin/bash
  
NCPU=10 #Number of CPUs

PARALLEL="/path/to/script/parallel.sh" #e.g. /data/tanoramb/Scripts/other/parallel.sh in Astro

RFMIX15="/path/to/Software/RFmix/RFMix_v1.5.4/PopPhased/RFMix_PopPhased"  #e.g. /data/tanoramb/Software/RFmix/RFMix_v1.5.4/PopPhased/RFMix_PopPhased in Astro
DIR="/path/to/data/rfmix15_files" #Path to directory containing files for RFMix 1.5 generated by Master script of previous step.

${PARALLEL} -j ${NCPU} -r "${RFMIX15} -a ${DIR}/rfmix15.chr*.rfmix15.hap -p ${DIR}/rfmix15.chr*.rfmix15.class -m ${DIR}/pos_chr*.cm.txt -fb 1 -w 0.2 -o rfmix15.out.chr* > rfmix15.chr*.log 2>&1" 1 22 2 21 3 20 4 19 5 18 6 17 7 16 8 15 9 14 10 13 11 12

#Gzipping heavy files
for CHR in $(seq 1 22); do
    gzip rfmix15.out.chr${CHR}.0.ForwardBackward.txt
    gzip rfmix15.out.chr${CHR}.allelesRephased0.txt
done
```


## 4. Further calculations

In the case of running RFMix 2.0, this section describes how to obtain the rephasing of original data and to get Viterbi files in the format of RFMix 1.5.

For the Re-phasing the suite [Tractor](https://github.com/Atkinson-Lab/Tractor) is used.

```
#! /bin/bash
  
 
DATA="/path/to/data" #Path to directory containing files and subdirectories generated by Master script in step 2
RFMIX2DIR="/path/to/rfmix2_data" #Path to directory containing files generated by RFMix 2.0 in step 3
UKMSP="python /path/to/Software/Tractor/UnkinkMSPfile.py" #e.g. /data/tanoramb/Software/Tractor/UnkinkMSPfile.py in Astro
UKGEN="python /path/to/Software/Tractor/UnkinkGenofile.py" #e.g. /data/tanoramb/Software/Tractor/UnkinkGenofile.py in Astro

for CHR in $(seq 1 22); do
    cp -s ${RFMIX2DIR}/rfmix20.chi.chr${CHR}.deconvoluted* .
    echo "Unkinking MSP file - CHR ${CHR}"
    ${UKMSP} --msp rfmix20.chi.chr${CHR}.deconvoluted
    echo "Gunzipping VCF CHR ${CHR}"
    gunzip -c ${DATA}/phased.intersected.norm.chi.chr${CHR}.vcf.gz > phased.intersected.norm.chi.chr${CHR}.vcf
    echo "Unkinking VCF CHR ${CHR}"
    ${UKGEN} --switches rfmix20.chi.chr${CHR}.deconvoluted.switches.txt --genofile phased.intersected.norm.chi.chr${CHR}
    echo "DONE CHR ${CHR}"
done
```

For getting Viterbi files in the RFMix 1.5 format from RFMix 2.0 output:
```
#! /bin/bash
  
DIR="/path/to/rfmix2_data" #Path to directory containing files generated by RFMix 2.0 in step 3

head -n 1 ${DIR}/rfmix20.chi.chr1.deconvoluted.msp.tsv | cut -d":" -f2 | awk '{print $1"\n"$2"\n"$3}' > rfmix20.Viterbi.coding.txt
sed -n 2p ${DIR}/rfmix20.chi.chr1.deconvoluted.msp.tsv | cut -f7- | sed 's/\t/\n/g' > rfmix20.Viterbi.columns.txt

for CHR in $(seq 1 22); do
    tail -n +3 ${DIR}/rfmix20.chi.chr${CHR}.deconvoluted.msp.tsv | cut -f7- | sed 's/\t/ /g' | sed 's/$/ /g' > rfmix20.out.chr${CHR}.0.Viterbi.txt
done
```

Any request at tanoramb@gmail.com


