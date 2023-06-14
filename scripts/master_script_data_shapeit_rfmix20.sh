#! /bin/bash

TSTART=$(date +%s)

CHR=0

CMD=$(echo "${0} ${*}")

function usage() {
    echo "Usage: $0 -c=<int>"
    echo "  -c,--chromosome         Number of the chromosome to generate data to. Should be an int 1 to 22."
    exit 1
}

if [ $# -eq 0 ]; then
    usage
else
    for i in "$@"; do
        case $i in
            -c=*|--chromosome=*) CHR="${i#*=}"; shift;;
            -p=*|--parameters=*) PARAMS="${i#*=}"; shift;;
            *)
                echo "Bad argument"
                usage
            ;;
        esac
    done
fi

if [[ "${CHR}" -lt 1 || "${CHR}" -gt 22 ]]; then
    echo "ERROR: Chromosome must be a number between 1 and 22"
    exit 1
fi

BEAGLE5="/data/tanoramb/Software/Beagle/beagle.12Jul19.0df.jar"
BEGLE5CMD="java -Xmx10g -jar ${BEAGLE5}"
SHAPEIT2="/data/tanoramb/Software/Shapeit/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
INTERSEC="/data/tanoramb/Scripts/v2/vcf_intersection_wrefsnps_resolvestrissues.v2.py"
TABIX="/data/tanoramb/Software/Tabix/tabix-0.2.6/tabix"
PLINK1="/data/tanoramb/Software/Plink1/plink"
PLINK2="/data/tanoramb/Software/Plink2/plink2"
VCFMERGE="/data/tanoramb/Software/VCFTools/vcftools-vcftools-15db94b/src/perl/vcf-merge"

MAPDIR="/data/tanoramb/Data/recombinationmaps/beagle5"
MAPDIR2="/data/tanoramb/Data/recombinationmaps/impute2"
SHAPEITDATA="/data/tanoramb/Data/shapeit_ref/1000GP_Phase3" #GOT FROM https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html

CHI="/data/tanoramb/Data/chi271"
PEL="/data/tanoramb/Data/peruvian"
IBS="/data/tanoramb/Data/iberian"
YRI="/data/tanoramb/Data/yoruba"

mkdir -p shapeit2_files

POPS=("PEL" "IBS" "YRI")

SETS=()

gunzip -c ${CHI}/chi271.chr${CHR}.rmd.vcf.gz | sed '/^#/d' | cut -f2 > onlypos_chr${CHR}.txt
gunzip -c ${IBS}/1KG_IBS_chr${CHR}.recode.rmd.vcf.gz | sed '/^#/d' | grep -w -f onlypos_chr${CHR}.txt | cut -f1,2,4,5 > pos_chr${CHR}.txt
cut -f2 pos_chr${CHR}.txt > .tmp.pos_chr${CHR}

gunzip -c ${CHI}/chi271.chr${CHR}.rmd.vcf.gz | sed -n '/^#/p' > chi.chr${CHR}.vcf
gunzip -c ${CHI}/chi271.chr${CHR}.rmd.vcf.gz | sed '/^#/d' | grep -w -f .tmp.pos_chr${CHR} >> chi.chr${CHR}.vcf
SETS+=("chi.chr${CHR}.vcf")

for POP in ${POPS[*]}; do
    name=$(echo ${POP} | awk '{print tolower($0)}')
    if [[ "${POP}" == "PEL" ]]; then
        gunzip -c ${!POP}/1KG_${POP}_chr${CHR}.recode.rmd.ge95.vcf.gz | sed -n '/^#/p' > ${name}.chr${CHR}.vcf
        gunzip -c ${!POP}/1KG_${POP}_chr${CHR}.recode.rmd.ge95.vcf.gz | grep -w -f .tmp.pos_chr${CHR} >> ${name}.chr${CHR}.vcf
    fi
    if [[ "${POP}" == "IBS" || "${POP}" == "YRI" ]]; then
        gunzip -c ${!POP}/1KG_${POP}_chr${CHR}.recode.rmd.vcf.gz | sed -n '/^#/p' > ${name}.chr${CHR}.vcf
        gunzip -c ${!POP}/1KG_${POP}_chr${CHR}.recode.rmd.vcf.gz | grep -w -f .tmp.pos_chr${CHR} >> ${name}.chr${CHR}.vcf
    fi
    SETS+=("${name}.chr${CHR}.vcf")
done


${INTERSEC} -s pos_chr${CHR}.txt -i $(echo ${SETS[@]} | tr ' ' ',') -o common_positions.chr${CHR}.txt -p > localancestry.data.chr${CHR}.log 2>&1

bgzip intersected.norm.chi.chr${CHR}.vcf
${TABIX} -p vcf intersected.norm.chi.chr${CHR}.vcf.gz

${SHAPEIT2} --input-vcf intersected.norm.chi.chr${CHR}.vcf.gz --input-map ${SHAPEITDATA}/genetic_map_chr${CHR}_combined_b37.txt --input-ref ${SHAPEITDATA}/1000GP_Phase3_chr${CHR}.hap.gz ${SHAPEITDATA}/1000GP_Phase3_chr${CHR}.legend.gz ${SHAPEITDATA}/1000GP_Phase3.sample -O shapeit2_files/phased.intersected.norm.chi.chr${CHR} >> localancestry.data.chr${CHR}.log 2>&1

${SHAPEIT2} -convert --input-haps shapeit2_files/phased.intersected.norm.chi.chr${CHR} --output-vcf phased.intersected.norm.chi.chr${CHR}.vcf >> localancestry.data.chr${CHR}.log 2>&1
bgzip phased.intersected.norm.chi.chr${CHR}.vcf
${TABIX} -p vcf phased.intersected.norm.chi.chr${CHR}.vcf.gz

ALL=("phased.intersected.norm.chi.chr${CHR}.vcf.gz")
REFS=()
for POP in ${POPS[*]}; do
    name=$(echo ${POP} | awk '{print tolower($0)}')
    bgzip intersected.norm.${name}.chr${CHR}.vcf
    mv intersected.norm.${name}.chr${CHR}.vcf.gz phased.intersected.norm.${name}.chr${CHR}.vcf.gz
    ${TABIX} -p vcf phased.intersected.norm.${name}.chr${CHR}.vcf.gz
    ALL+=("phased.intersected.norm.${name}.chr${CHR}.vcf.gz")
    REFS+=("phased.intersected.norm.${name}.chr${CHR}.vcf.gz")
done

${VCFMERGE} -d -t -s -c none -R "0|0" ${ALL[@]} > phased.intersected.norm.ALL.chr${CHR}.vcf
bgzip phased.intersected.norm.ALL.chr${CHR}.vcf
${TABIX} -p vcf phased.intersected.norm.ALL.chr${CHR}.vcf.gz
${VCFMERGE} -d -t -s -c none -R "0|0" ${REFS[@]} > phased.intersected.norm.REFS.chr${CHR}.vcf
bgzip phased.intersected.norm.REFS.chr${CHR}.vcf
${TABIX} -p vcf phased.intersected.norm.REFS.chr${CHR}.vcf.gz


TEND=$(date +%s)
RUNTIME=$((TEND-TSTART))
echo "DONE: ${RUNTIME} sec" >> localancestry.data.chr${CHR}.log 2>&1
