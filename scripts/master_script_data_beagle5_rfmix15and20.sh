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
#BEGLE5CMD="java -Xmx10g -jar ${BEAGLE5}"
BEGLE5CMD="java -Xmx15g -jar ${BEAGLE5}"
SHAPEIT2="/data/tanoramb/Software/Shapeit/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
INTERSEC="/data/tanoramb/Scripts/v2/vcf_intersection_wrefsnps_resolvestrissues.v3.py"
TABIX="/data/tanoramb/Software/Tabix/tabix-0.2.6/tabix"
FORMATRFMIX15="/data/tanoramb/Scripts/v2/format_vcf_to_rfmix15.py"
PLINK1="/data/tanoramb/Software/Plink1/plink"
PLINK2="/data/tanoramb/Software/Plink2/plink2"
LAMPLPGEN="/data/tanoramb/Scripts/v2/phased_vcf_to_lampld_genfile.py"
LAMPLPHAP="/data/tanoramb/Scripts/v2/phased_vcf_to_lampld_hapfile.py"
VCFMERGE="/data/tanoramb/Software/VCFTools/vcftools-vcftools-15db94b/src/perl/vcf-merge"

MAPDIR="/data/tanoramb/Data/recombinationmaps/beagle5"
MAPDIR2="/data/tanoramb/Data/recombinationmaps/impute2"
#GOT FROM https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html
SHAPEITDATA="/data/tanoramb/Data/shapeit_ref/1000GP_Phase3"

CHI="/data/tanoramb/Data/chi271"
PEL="/data/tanoramb/Data/peruvian"
IBS="/data/tanoramb/Data/iberian"
YRI="/data/tanoramb/Data/yoruba"

mkdir -p rfmix15_files

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
        gunzip -c ${!POP}/1KG_${POP}_chr${CHR}.recode.rmd.ge90.vcf.gz | sed -n '/^#/p' > ${name}.chr${CHR}.vcf
        gunzip -c ${!POP}/1KG_${POP}_chr${CHR}.recode.rmd.ge90.vcf.gz | grep -w -f .tmp.pos_chr${CHR} >> ${name}.chr${CHR}.vcf
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

${BEGLE5CMD} map=${MAPDIR}/plink.chr${CHR}.GRCh37.map gt=intersected.norm.chi.chr${CHR}.vcf.gz out=phased.intersected.norm.chi.chr${CHR} >> localancestry.data.chr${CHR}.log 2>&1
gunzip phased.intersected.norm.chi.chr${CHR}.vcf.gz
bgzip phased.intersected.norm.chi.chr${CHR}.vcf
${TABIX} -p vcf phased.intersected.norm.chi.chr${CHR}.vcf.gz
gunzip -k phased.intersected.norm.chi.chr${CHR}.vcf.gz

${FORMATRFMIX15} -i phased.intersected.norm.chi.chr${CHR}.vcf -o phased.intersected.norm.chi.chr${CHR}.rfmix15.hap >> localancestry.data.chr${CHR}.log 2>&1

CLASS=0
ln -s phased.intersected.norm.chi.chr${CHR}.rfmix15.hap .tmp.chi.${CLASS}.${CHR}.hap
grep "CHROM" phased.intersected.norm.chi.chr${CHR}.vcf | cut -f10- | tr '\t' '\n' | awk -v c=${CLASS} 'BEGIN{s=""}{s=s c " " c " "}END{print substr(s,1,length(s)-1)}' > phased.intersected.norm.chi.chr${CHR}.rfmix15.class
ln -s phased.intersected.norm.chi.chr${CHR}.rfmix15.class .tmp.chi.${CLASS}.${CHR}.class
${PLINK2} --vcf phased.intersected.norm.chi.chr${CHR}.vcf.gz --make-bed --out phased.intersected.norm.chi.chr${CHR} >> localancestry.data.chr${CHR}.log 2>&1
${PLINK1} --bfile phased.intersected.norm.chi.chr${CHR} --cm-map ${MAPDIR2}/chr@.impute2.map --make-bed --out phased.intersected.norm.chi.chr${CHR}.cm >> localancestry.data.chr${CHR}.log 2>&1
cut -f3 phased.intersected.norm.chi.chr${CHR}.cm.bim > rfmix15_files/pos_chr${CHR}.cm.txt
\rm -rf phased.intersected.norm.chi.chr${CHR}.vcf
HAPSET+=(".tmp.chi.${CLASS}.${CHR}.hap")
CLSSET+=(".tmp.chi.${CLASS}.${CHR}.class")

ALL=("phased.intersected.norm.chi.chr${CHR}.vcf.gz")
REFS=()
CLASS=1
for POP in ${POPS[*]}; do
    name=$(echo ${POP} | awk '{print tolower($0)}')
    bgzip intersected.norm.${name}.chr${CHR}.vcf
    ## 1000G genomes are already phased so Beagle.v5 call is skipped
    cp intersected.norm.${name}.chr${CHR}.vcf.gz phased.intersected.norm.${name}.chr${CHR}.vcf.gz
    #${TABIX} -p vcf intersected.norm.${name}.chr${CHR}.vcf.gz
    #${BEGLE5CMD} map=${MAPDIR}/plink.chr${CHR}.GRCh37.map gt=intersected.norm.${name}.chr${CHR}.vcf.gz out=phased.intersected.norm.${name}.chr${CHR} >> localancestry.data.chr${CHR}.log 2>&1
    #gunzip phased.intersected.norm.${name}.chr${CHR}.vcf.gz
    #bgzip phased.intersected.norm.${name}.chr${CHR}.vcf
    ${TABIX} -p vcf phased.intersected.norm.${name}.chr${CHR}.vcf.gz
    gunzip -k phased.intersected.norm.${name}.chr${CHR}.vcf.gz
    ${FORMATRFMIX15} -i phased.intersected.norm.${name}.chr${CHR}.vcf -o phased.intersected.norm.${name}.chr${CHR}.rfmix15.hap >> localancestry.data.chr${CHR}.log 2>&1
    ln -s phased.intersected.norm.${name}.chr${CHR}.rfmix15.hap .tmp.${name}.${CLASS}.${CHR}.hap
    grep "CHROM" phased.intersected.norm.${name}.chr${CHR}.vcf | cut -f10- | tr '\t' '\n' | awk -v c=${CLASS} 'BEGIN{s=""}{s=s c " " c " "}END{print substr(s,1,length(s)-1)}' > phased.intersected.norm.${name}.chr${CHR}.rfmix15.class
    ln -s phased.intersected.norm.${name}.chr${CHR}.rfmix15.class .tmp.${name}.${CLASS}.${CHR}.class
    \rm -rf phased.intersected.norm.${name}.chr${CHR}.vcf
    HAPSET+=(".tmp.${name}.${CLASS}.${CHR}.hap")
    CLSSET+=(".tmp.${name}.${CLASS}.${CHR}.class")

    ALL+=("phased.intersected.norm.${name}.chr${CHR}.vcf.gz")
    REFS+=("phased.intersected.norm.${name}.chr${CHR}.vcf.gz")

    ((CLASS+=1))
done
paste -d " " ${CLSSET[@]} > rfmix15_files/rfmix15.chr${CHR}.rfmix15.class
paste -d "" ${HAPSET[@]} > rfmix15_files/rfmix15.chr${CHR}.rfmix15.hap


${VCFMERGE} -d -t -s -c none -R "0|0" ${ALL[@]} > phased.intersected.norm.ALL.chr${CHR}.vcf
bgzip phased.intersected.norm.ALL.chr${CHR}.vcf
${TABIX} -p vcf phased.intersected.norm.ALL.chr${CHR}.vcf.gz
${VCFMERGE} -d -t -s -c none -R "0|0" ${REFS[@]} > phased.intersected.norm.REFS.chr${CHR}.vcf
bgzip phased.intersected.norm.REFS.chr${CHR}.vcf
${TABIX} -p vcf phased.intersected.norm.REFS.chr${CHR}.vcf.gz


TEND=$(date +%s)
RUNTIME=$((TEND-TSTART))
echo "DONE: ${RUNTIME} sec" >> localancestry.data.chr${CHR}.log 2>&1
