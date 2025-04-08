#!/bin/bash
# This script applies the GATK filters we use.

INPUT_DIR=$1
INTERVAL=$2
[[ -z $INTERVAL ]] && INTERVAL=${LSB_JOBINDEX}
INPUT_VCF=$(find ${INPUT_DIR} -name "*.${INTERVAL}.setmiss.annAB.vcf.gz")
COUNT_NAME="${INPUT_DIR##*/}.${INTERVAL}"

NAMEBASE="${INPUT_VCF%.*.*}"
COUNT_FILE="${NAMEBASE}.counts.csv"

# Set up counts file
echo "Name,Original,QUAL,LCR,no_mono,maxDP,ABHet,splitMA,geno" > "${COUNT_FILE}"
echo -n "${COUNT_NAME},$(zcat ${INPUT_VCF} | grep -v '^#' | wc -l)," >> "${COUNT_FILE}"

# Filter on GLNexus QUAL score (AQ tag)
bcftools view \
	-e "INFO/AQ < 100" \
	-Oz -o ${NAMEBASE}.QUAL.vcf.gz \
	${INPUT_VCF}

echo -n "$(zcat ${NAMEBASE}.QUAL.vcf.gz | grep -v '^#' | wc -l)," >> "${COUNT_FILE}"

# Filter LCR
LCR="/scripts/LCR-hs38.bed"
bedtools subtract \
	-header \
	-A -a ${NAMEBASE}.QUAL.vcf.gz \
	-b ${LCR} \
	| bgzip -c > ${NAMEBASE}.QUAL.LCR.vcf.gz 

tabix ${NAMEBASE}.QUAL.LCR.vcf.gz \
&& rm ${NAMEBASE}.QUAL.vcf.gz*

echo -n "$(zcat ${NAMEBASE}.QUAL.LCR.vcf.gz | grep -v '^#' | wc -l)," >> "${COUNT_FILE}"

# Filter monomorphics
bcftools view \
    -e ' INFO/AF = 0 || INFO/AF = 1' \
    -Oz -o ${NAMEBASE}.QUAL.LCR.nomono.vcf.gz \
    ${NAMEBASE}.QUAL.LCR.vcf.gz \
&& rm ${NAMEBASE}.QUAL.LCR.vcf.gz*

tabix ${NAMEBASE}.QUAL.LCR.nomono.vcf.gz

echo -n "$(zcat ${NAMEBASE}.QUAL.LCR.nomono.vcf.gz | grep -v '^#' | wc -l)," >> "${COUNT_FILE}"

# Filter variants with mean read depth over 500
vcftools \
	--gzvcf ${NAMEBASE}.QUAL.LCR.nomono.vcf.gz \
	--max-meanDP 500 \
	--recode-INFO-all \
	--recode \
	--stdout \
	| bgzip -c > ${NAMEBASE}.QUAL.LCR.nomono.maxDP.vcf.gz
[[ -f ${NAMEBASE}.QUAL.LCR.nomono.maxDP.vcf.gz ]] \
|| bcftools head ${NAMEBASE}.QUAL.LCR.nomono.vcf.gz \
	| bgzip -c > ${NAMEBASE}.QUAL.LCR.nomono.maxDP.vcf.gz

tabix ${NAMEBASE}.QUAL.LCR.nomono.maxDP.vcf.gz

[[ -f ${NAMEBASE}.QUAL.LCR.nomono.maxDP.vcf.gz ]] && rm ${NAMEBASE}.QUAL.LCR.nomono.vcf.gz*

echo -n "$(zcat ${NAMEBASE}.QUAL.LCR.nomono.maxDP.vcf.gz | grep -v '^#' | wc -l)," >> "${COUNT_FILE}"

# Filter variants outside 0.25 < ABHet < 0.75
${GATK} \
    --java-options "-Xmx8g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	SelectVariants \
	-R ${REF_FASTA} \
	-V ${NAMEBASE}.QUAL.LCR.nomono.maxDP.vcf.gz \
	-select '(!vc.hasAttribute("ABHet")) || (ABHet >= 0.25 && ABHet <= 0.75)' \
	-O ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.vcf.gz

[[ -f ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.vcf.gz ]] && rm ${NAMEBASE}.QUAL.LCR.nomono.maxDP.vcf.gz*

echo -n "$(zcat ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.vcf.gz | grep -v '^#' | wc -l)," >> "${COUNT_FILE}"

bcftools norm \
	-m -any \
	${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.vcf.gz \
| bcftools annotate \
	--set-id '%CHROM\:%POS\:%REF\:%ALT' \
| bgzip -c > ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.vcf.gz
tabix -f ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.vcf.gz

[[ -f ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.vcf.gz ]] && rm ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.vcf.gz*

echo -n "$(zcat ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.vcf.gz | grep -v '^#' | wc -l)," >> "${COUNT_FILE}"

# Apply genotyping rate filter
plink2 --vcf ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.vcf.gz \
	--double-id \
	--geno 0.05 \
	--vcf-half-call missing \
	--memory 9000 \
	--threads 1 \
	--make-just-bim \
	--out ${NAMEBASE}.geno

if [[ -s ${NAMEBASE}.geno.bim ]]; then
# Get a list of positions to keep
cat ${NAMEBASE}.geno.bim \
	| cut -f2 \
	| awk -F ":" '{print $1,$2,$2}' OFS="\t" \
	> ${NAMEBASE}.geno.keep_regions.bed

# Filter variants
bedtools intersect \
	-a ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.vcf.gz \
	-b ${NAMEBASE}.geno.keep_regions.bed \
	-header \
	| bgzip -c > ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.geno.vcf.gz
	tabix -f ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.geno.vcf.gz
	SUCCESS=$?
else 
	mv ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.vcf.gz ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.geno.vcf.gz
	mv ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.vcf.gz.tbi ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.geno.vcf.gz.tbi
	SUCCESS=$?
fi

echo -n "$(zcat ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.geno.vcf.gz | grep -v '^#' | wc -l)" >> "${COUNT_FILE}"

[[ $SUCCESS -eq 0 ]] \
&& rm ${NAMEBASE}.QUAL.LCR.nomono.maxDP.ABHet.splitMA.vcf.gz* \
&& rm ${INPUT_VCF}* \
&& rm ${NAMEBASE}.geno.{log,bim,keep_regions.bed}

exit $SUCCESS
