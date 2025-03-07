#!/bin/bash
RECAL_VCF="$1"
NAMEBASE="${RECAL_VCF%.*.*}"
COUNT_FILE="${NAMEBASE}.counts.csv"

echo "Recal,${NAMEBASE}.vcf.gz,$(zcat ${NAMEBASE}.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE

${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	SelectVariants \
	-R ${REF_FASTA} \
	-L ${CHR} \
	-V ${NAMEBASE}.vcf.gz \
	--exclude-filtered \
	-O ${NAMEBASE}-PASS.vcf.gz

echo "Filtered,${NAMEBASE}-PASS.vcf.gz,$(zcat ${NAMEBASE}-PASS.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE
	
LCR="/scripts/LCR-hs38.bed"
bedtools subtract -header -A -a ${NAMEBASE}-PASS.vcf.gz -b ${LCR} | bgzip -c > ${NAMEBASE}-PASS-LCR.vcf.gz 
bcftools index -t --threads $LSB_MAX_NUM_PROCESSORS ${NAMEBASE}-PASS-LCR.vcf.gz

echo "LCR Filtered,${NAMEBASE}-PASS-LCR.vcf.gz,$(zcat ${NAMEBASE}-PASS-LCR.vcf.gz | grep -v '#' | wc -l)" >> $COUNT_FILE

${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	SelectVariants \
	-R ${REF_FASTA} \
	-L ${CHR} \
	-V ${NAMEBASE}-PASS-LCR.vcf.gz \
	--exclude-non-variants \
	-O ${NAMEBASE}-PASS-LCR-nonVariants.vcf.gz

echo "Non-Variant Filtered,${NAMEBASE}-PASS-LCR-nonVariants.vcf.gz,$(zcat ${NAMEBASE}-PASS-LCR-nonVariants.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE

${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	SelectVariants \
	-R ${REF_FASTA} \
	-L ${CHR} \
	-V ${NAMEBASE}-PASS-LCR-nonVariants.vcf.gz \
	-select 'vc.hasGenotypes() && vc.getCalledChrCount(vc.getAltAlleleWithHighestAlleleCount())/(1.0*vc.getCalledChrCount()) < 1.0'  \
	-O ${NAMEBASE}-PASS-LCR-nonVariants-AF1.vcf.gz

echo "AF Filtered,${NAMEBASE}-PASS-LCR-nonVariants-AF1.vcf.gz,$(zcat ${NAMEBASE}-PASS-LCR-nonVariants-AF1.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE

java -Xmx80g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -XX:+UseSerialGC \
	-jar /ref/GATK360.jar \
	-T VariantAnnotator \
	-A AlleleBalance \
	-R /ref/20190812_GATK_38_googlebundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
	-V ${NAMEBASE}-PASS-LCR-nonVariants-AF1.vcf.gz \
	-o ${NAMEBASE}-PASS-LCR-nonVariants-AF1-ABannotated.vcf.gz

echo "AB Annotated,${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1-ABannotated.vcf.gz,$(zcat ${NAMEBASE}-PASS-LCR-nonVariants-AF1-ABannotated.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE

${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	SelectVariants \
	-R ${REF_FASTA} \
	-L ${CHR} \
	-V ${NAMEBASE}-PASS-LCR-nonVariants-AF1-ABannotated.vcf.gz \
	-select '(!vc.hasAttribute("ABHet")) || (ABHet >= 0.30 && ABHet <= 0.70)' \
	-O ${NAMEBASE}-PASS-LCR-nonVariants-AF1-ABfiltered.vcf.gz
SUCCESS=$?

echo "AB Filtered,${NAMEBASE}-PASS-LCR-nonVariants-AF1-ABfiltered.vcf.gz,$(zcat ${NAMEBASE}-PASS-LCR-nonVariants-AF1-ABfiltered.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE
exit $SUCCESS